# Main function ------------------------------------------------------------------------
globalVariables(c("gr","colts","cgr","cid","cms","ms","kmeans"))

#' Clustering Algorithm to HDLSS data.
#'
#' @param data A numeric matrix or data frame with all numeric columns. If a matrix or data frame, rows correspond to variables (d) and columns correspond to observations (n).
#' @param alpha A real number in the range (0, 1) indicanting the threshold parameter to be compared with p-values in the clustering procedure.
#' @param n.cores A number processor cores (see detectCores).
#' @param ... not used.
#' 
#' @return Results
#' 
#' A vector of integers indicating the cluster to which each variable is allocated and a cluster cross table of P-values.
#' 
#' @examples 
#' #data(diflogadenoma) # loads data
#' #cl <- ppclust_h(diflogadenoma[, 10:13], n.cores = 2)
#' #cl
#' @importFrom stats aggregate cov median pnorm reshape var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel clusterExport makeCluster parApply stopCluster
#' @export

#criar funcao para otimizar o agrup/tomar o menor numero de grupos

ppclust_h <- function(data,
                      alpha,
                      n.cores = 1,
                      ...)
{
  
  
  anovaRank <- function(x)
  {
    sigma4 <- x$sigma4
    a <- nrow(x)
    Ri. <- x$mean
    R.. <- mean(Ri.)
    MSTr <- 1 / (a - 1) * sum((Ri. - R..) ** 2)
    S2ri <- x$v22p
    MSEr <- 1 / a * sum(S2ri)
    Fr <- MSTr / MSEr
    est <- sqrt(a) * (Fr - 1)
    v22 <- 1 / a * sum(S2ri)
    tal2p <- x$tal2p
    tal2 <- 1 / a * sum(tal2p)
    var <- tal2 / v22 ** 2
    pv <- 1 - pnorm(est / sqrt(var))
    return(pv)
  }
  
  sortMean <- function(x)
  {
    x[order(x$mean), ]
  }
  
  computeSigma4 <- function(x)
  {
    n <- length(x)
    if (n > 1) {
      s4hat <- (1 / n * sum((x - mean(x)) ** 2)) ** 2
      s4jack <- n * s4hat
      m <- matrix(c(rep(x, n)), byrow = T, ncol = n)
      diag(m) <- NA
      s4hatni <- apply(m, 1, function(x) {
        xi = x[!is.na(x)]
        (sum((xi - mean(xi)) ** 2 / (n - 1))) ** 2
      })
      s4hatn <- sum(s4hatni)
      sigma4 <- s4jack - (n - 1) / n * s4hatn
    } else{
      sigma4 <- 0
    }
    sigma4
  }
  
  partition <- function(data, cts, cms, cl.cores)
  {
    if (nrow(data) < 2) {
      return(list(pt = 1, cts = cts))
    } else{
      datax <- data$mean
      #iter.m <- floor(sqrt(length(datax)))
      partitions <- suppressWarnings(kmeans(datax, length(cts), algorithm = "Hartigan-Wong"))
      return(list(pt = partitions$cluster, cts = c(partitions$centers)))
    }
  }
  
  groupAnovaRank <- function(data, alpha)
  {
    
    p_values <- function(data){
      n <- length(unique(data$gr))
      m <- matrix(NA, nrow = n, ncol = n)
      indexes <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
      pvs <- apply(indexes, 1, function(x) anovaRank(data[data$gr == x[1] | data$gr == x[2], ]))
      for(i in 1:nrow(indexes)) m[indexes[i,1], indexes[i,2]] <- pvs[i]
      return(m)
    }
    
    m <- p_values(data)
    mpv <- max(m, na.rm = T)
    while (mpv >= alpha) {
      idxs <- c(which(m == max(m, na.rm = T), arr.ind = T))
      x <- data$gr
      n <- length(unique(data$gr))
      x[x == max(idxs)] <- n + 1
      x[x == min(idxs)] <- n + 1
      x <- as.numeric(factor(x))
      data$gr <- x
      m <- p_values(data)
      mpv <- max(m, na.rm = T)
    }
    
    n <- length(unique(data$gr))
    m <- p_values(data)
    indexes <- cbind(1:n, 1:n)
    pvs <- apply(indexes, 1, function(x) anovaRank(data[data$gr == x[1] | data$gr == x[2], ]))
    diag(m) <- pvs
    m <- as.data.frame(matrix(format.pval(c(m), digits = 6), nrow = n, ncol = n))
    colnames(m) <- 1:n
    list(cluster = data$gr, pv.matrix = m)
  }
  
  fib <- function(fibseq, x)
  {
    c(fibseq, fibseq[x - 1] + fibseq[x - 2])
  }
  
  if(!is.data.frame(data) & !is.matrix(data))
    stop('the data should be a numeric matrix or data frame with all numeric columns')
  
  if (any(apply(data, 2, is.numeric) == FALSE))
    stop('the data contains non-numeric values')
  
  if (!(alpha > 0 & alpha < 1))
    stop("no valid value for 'alpha'. For help see ?ppclust")
  
  sample.size <- apply(data, 1, function(x) length(x[!is.na(x)]))
  
  if(any(sample.size < 2))
    stop("all sample sizes must be larger than 2")
  
  cl.cores <- parallel::makeCluster(n.cores)
  rankData <- matrix(rank(c(as.matrix(data)), na.last = 'keep'),
                     nrow = nrow(data), ncol = ncol(data))
  ncole <- ncol(rankData)
  nrowe <- nrow(rankData)
  environment(computeSigma4) <- .GlobalEnv
  parallel::clusterExport(cl.cores, varlist = c("computeSigma4"), envir = environment())
  
  X <- parallel::parApply(cl.cores, rankData, 1, function(x) {
    x1 <- x[!is.na(x)]
    c(computeSigma4(x1), var(x1), length(x1), mean(x1))
  })
  X <- t(X)
  sigma4 <- X[, 1]
  S2ri <- X[, 2]
  n <- X[, 3]
  means <- X[, 4]
  tal2p <- 2 / (n * (n - 1)) * sigma4
  v22p <- 1 / n * S2ri
  
  ppData <- data.frame(id = 1:nrowe,
                       mean = means,
                       sigma4 = sigma4,
                       v22p = v22p,
                       tal2p = tal2p,
                       gr = 99999)
  cms <- 2
  cid <- 1
  cs4 <- 3
  cv22 <- 4
  ctal2 <- 5
  cgr <- 6
  
  fstart <- 1
  fend <- nrowe
  g <- 1
  
  pValue <- anovaRank(ppData)
  
  if (pValue > alpha) {
    parallel::stopCluster(cl.cores)
    ppData$gr <- g
    return(list(cluster = ppData$gr,
                P.values = pValue))
  } else{
    
    ppData <- sortMean(ppData)
    mid <- ppData$mean[c(1, fend)]
    new.mid <- c()
    k <- 3
    fibseq <- c(1, 1)
    for (i in 1:15) fibseq <- fib(fibseq, i)
    
    cat('clustering ... \n')
    pb <- txtProgressBar(min = 0, max = fend, style = 3)
    while (fstart <= fend) {
      mid <- matrix(c(mid, new.mid), ncol = 1)
      groups <- partition(ppData[fstart:fend, ], mid, cms, cl.cores)
      ptt <- groups$pt
      mid <- groups$cts
      mg <- 0
      pValue <- numeric()
      
      for (i in 1:length(unique(ptt))) {
        pValue[i] <- anovaRank(ppData[c(fstart - 1 + which(ptt == i)), ])
        if (pValue[i] > alpha | is.nan(pValue[i])) {
          ppData$gr[c(fstart - 1 + which(ptt == i))] <- g
          g <- g + 1
          ppData <- ppData[order(ppData$gr), ]
          mg <- 1
        }
      }
      fstart <- sum(ppData$gr != 99999) + 1
      if (mg != 0) {
        mid <- c()
        k <- 3
      } else{
        mid <- c()
        k <- k + 1
      }
      if (fstart < fend)
        new.mid <- ppData[fstart:(fstart - 1 + fibseq[k]), "mean"]
    }
  }
  
  parallel::stopCluster(cl.cores)
  resul <- groupAnovaRank(ppData, alpha)
  setTxtProgressBar(pb, fend)
  close(pb)
  
  ppData$gr <- resul$cluster
  
  ppData <- ppData[order(ppData$gr), ]
  groupclass <- ppData[, c("id", "mean", "gr")]
  groupclass0 <- groupclass[which(groupclass$gr == 0), ]
  groupclass <- groupclass[which(groupclass$gr != 0), ]
  groupc <- tapply(groupclass$mean, groupclass$gr, function(x) mean(x))
  groupc <- as.numeric(names(sort(groupc)))
  leng <- length(groupc)
  for(i in 1:leng){
    groupclass[groupclass$gr == groupc[i], ] <- leng + i
  }
  
  ppData$gr <- c(groupclass0$gr, as.numeric(as.factor(groupclass$gr)))
  ppData <- ppData[order(ppData$id), ]
  
  return(list(cluster = ppData$gr,
              P.values = resul$pv.matrix))
}

