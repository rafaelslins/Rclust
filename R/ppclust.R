
# Main function ------------------------------------------------------------------------
globalVariables(c("data.table","as.data.table"))

#' @title Clustering Algorithm to HDLSS data.
#' @description Clustering Algorithm to HDLSS data.

#' @usage ppclust(data, alpha, ...) 
#' @param data A numeric matrix or data frame with all numeric columns. If a matrix or data frame, rows correspond to variables (d) and columns correspond to observations (n).
#' @param dataset A numeric matrix or data frame with all numeric columns. If a matrix or data frame, rows correspond to variables (d) and columns correspond to observations (n).
#' @param alpha A real number in the range (0, 1) indicanting the threshold parameter to be compared with p-values in the clustering procedure.
#' @param ... not used.
#' 
#' @return Results
#' 
#' A vector of integers indicating the cluster to which each variable is allocated and a cluster cross table of P-values.
#' 
#' @examples 
#' #data(diflogadenoma) # loads data
#' #cl <- ppclust(diflogadenoma[, 10:13])
#' #cl
#' @importFrom stats aggregate cov median pnorm reshape var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

ppclust <- function(data,
                    alpha,
                    ...)
{
  
  ppclust <- function(dataset, alpha, ...) {
    
    if (any(apply(dataset, 2, is.numeric) == FALSE))
      stop('dataset contains non-numeric values')
    
    if (!(alpha > 0 & alpha < 1))
      stop("enter a valid 'alpha'. See ?ppclust.")
    
    if(!is.data.frame(data) & !is.matrix(data))
      stop('the data should be a numeric matrix or data frame with all numeric columns')
    
    if (any(apply(data, 2, is.numeric) == FALSE))
      stop('the data contains non-numeric values')
    
    if (!(alpha > 0 & alpha < 1))
      stop("no valid value for 'alpha'. For help see ?ppclust")
    
    sample.size <- apply(data, 1, function(x) length(x[!is.na(x)]))
    
    if(any(sample.size < 2))
      stop("all sample sizes must be larger than 2")
    
    rankData <- matrix(rank(c(as.matrix(data)), na.last = 'keep'),
                       nrow = nrow(data), ncol = ncol(data))
    rankData <- as.data.table(rankData)
    
    nrowe <- nrow(rankData)
    tstart <- 1
    tend <- nrowe
    g <- 1
    testpp <- 0
    means <- rowMeans(rankData, na.rm = T)
    vars <- apply(rankData, 1, function(x) var(x, na.rm = T))
    sigma4 <- computeSigma4(rankData)
    ppData <- data.table(id = 1:nrowe,
                         gr = 99999,
                         mean = means,
                         var = vars,
                         n = sample.size,
                         sigma4,
                         colts = 0)
    ppData <- sortMedian(rankData, ppData)
    ppData <- sortCenter(ppData)
    pValue <- anovaRank(ppData[tstart:tend, ])
    
    if (pValue > alpha) {
      ppData[, "gr"] <- 1
      return(list(cluster = ppData[, gr],
                  P.values = pValue))
    } else{
      
      cat('clustering ... \n')
      L <- (tend - tstart + 1) / 2
      tend <- tstart + floor(L - 1)
      pb <- txtProgressBar(min = 0, max = nrowe, style = 3)
      while (tstart <= nrowe) {
        if (tstart == nrowe) {
          ppData[tstart, "gr"] <- 0
          tstart <- tstart + 1
        } else{
          if (pValue > alpha) g <- g + 1
          
          ppData[, "colts"] <- 0
          j1 <- ppData[, colts]
          j1[which(tstart <= c(1:nrowe) & c(1:nrowe) <= tend)] <- 1
          ppData[, "colts"] <- j1
          
          n <- sum(ppData[, colts])
          
          if (n == 1){
            ppData[tstart, "gr"] <- 0
            tstart <- tend + 1
            tend <- nrowe
            testpp <- -1
          }
          
          ppDataSubset <- ppData[ppData[, colts] != 0, ]
          if (testpp == 0) pValue <- anovaRank(ppDataSubset)
          
          if (pValue > alpha) {
            j2 <- ppData[, gr]
            j2[which(tstart <= c(1:nrowe) & c(1:nrowe) <= tend)] <- g
            ppData[, "gr"] <- j2
            
            tstart <- tend + 1
            tend <- nrowe
            setTxtProgressBar(pb, tstart)
            resul <- indTest(ppData, tstart, nrowe, g, pValue, alpha)
            
            tstart <- resul$st
            ppData <- resul$data
            setTxtProgressBar(pb, tstart)
          } else{
            L <- (tend - tstart + 1) * 0.9
            tend <- tstart + floor(L - 1)
            testpp <- 0
          }
        }
      }
      setTxtProgressBar(pb, nrowe)
      close(pb)
    }
    
    ppData <- ppData[order(ppData$gr), ]
    groupclass <- ppData[, 1:3]
    groupclass0 <- groupclass[groupclass$gr == 0, ]
    groupclass <- groupclass[groupclass$gr != 0, ]
    groupc <- tapply(groupclass$mean, groupclass$gr, function(x) mean(x))
    groupc <- as.numeric(names(sort(groupc)))
    leng <- length(groupc)
    for(i in 1:leng){
      groupclass[groupclass$gr == groupc[i], ] <- leng + i
    }
    
    ppData[, "gr"] <- c(groupclass0$gr, as.numeric(as.factor(groupclass$gr)))
    ppData <- ppData[order(ppData$id), ]
    
    n <- length(unique(ppData$gr))
    st <- min(unique(ppData$gr))
    en <- max(unique(ppData$gr))
    m <- matrix(NA, nrow = n, ncol = n)
    indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1))) + (st - 1)
    pvs <- apply(indices, 1, function(x)
      anovaRank(ppData[ppData$gr == x[1] | ppData$gr == x[2],]))
    m[lower.tri(m)] <- pvs
    indices <- cbind(st:en, st:en)
    pvs <- apply(indices, 1, function(x)
      anovaRank(ppData[ppData$gr == x[1] | ppData$gr == x[2],]))
    diag(m) <- pvs
    m <- as.data.frame(matrix(format.pval(c(m), digits = 6), nrow = n, ncol = n))
    colnames(m) <- st:en
    row.names(m) <- colnames(m)
    return(list(cluster = ppData[, gr],
                P.values = m))
  }
  
  
  # Internal functions ------------------------------------------------------------------------
  
  anovaRank <- function(x)
  {
    meang <- x[, "mean"]
    varg <- x[, "var"]
    sigma4 <- x[, "sigma4"]
    nfac <- nrow(x)
    ng <- x[, "n"]
    meant <- sum(meang) / nfac
    mst <- sum((meang - meant) ** 2) / (nfac - 1)
    mse2 <- sum(varg / ng) / nfac
    fr2 <- mst / mse2
    tau2 <- sum(sigma4 / (ng * (ng - 1))) * 2 / nfac
    asyvar <- tau2 / mse2 ** 2
    pv <- 1 - pnorm(sqrt(nfac) * (fr2 - 1) / sqrt(asyvar))
    return(pv)
  }
  
  sortMedian <- function(rank.data, pp.data)
  {
    rank.data[is.na(rank.data) == T] <- nrow(rank.data) + 1000
    medians <- apply(rank.data, 1, function(x) median(x, na.rm = T))
    pp.data <- data.table(pp.data, median = medians)
    pp.data[order(medians), ]
  }
  
  sortCenter <- function(x)
  {
    nr <- nrow(x)
    j1 <- numeric(nr)
    c1 <- round(0.40 * nr, 0)
    c2 <- round(0.60 * nr, 0)
    j1[c(1:nr) < c1] <- 1
    j1[c(1:nr) > c2] <- 1
    x[order(j1, x[, "median"]), ]
  }
  
  computeSigma4 <- function(data)
  {
    sig4hat <- function(x, n) (sum((x - mean(x)) ** 2) / n) ** 2
    
    est <- function(x) {
      xi <- x[!is.na(x)]
      n <- length(xi)
      if (n > 1) {
        s4jack <- n * sig4hat(xi, n)
        for (i in 1:n) {
          xni <- xi[-i]
          s4hatni <- sig4hat(xni, (n - 1))
          s4jack <- s4jack - (n - 1) / n * s4hatni
        }
        sigma4 <- s4jack
      } else{
        sigma4 <- 0
      }
      return(sigma4)
    }
    sigma4 <- apply(data, 1, est)
    return(sigma4)
  }
  
  indTest <- function(data, st, nrowe, g, pValue, alpha)
  {
    data[, "colts"] <- 0
    j1 <- data[, colts]
    j1[data[, gr] == g] <- 1
    data[, "colts"] <- j1
    count <- 0
    
    if (st > nrowe)
      return(list(st = st, data = data))
    
    for (i in st:nrowe) {
      data[i, "colts"] <- 1
      n <- sum(data[, colts])
      data2 <- data[data[, colts] != 0, ]
      pValue <- anovaRank(data2)
      if (pvalue > alpha)
        data[i, "gr"] <- g
      else
        data[i, "colts"] <- 0
    }
    data <- data[order(data[, gr]), ]
    count <- sum(data[, gr] <= g)
    st <- count + 1
    list(st = st, data = as.data.table(data))
  }
  
  