
# Main function ------------------------------------------------------------------------
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

ppclust_h <- function(data,
                      alpha,
                      n.cores = 1,
                      ...)
{
  
  if(!is.data.frame(data) & !is.matrix(data))
    stop('the data should be a numeric matrix or data frame with all numeric columns')
  
  if (any(apply(data, 2, is.numeric) == FALSE))
    stop('the data contains non-numeric values')
  
  if (!(alpha > 0 & alpha < 1))
    stop("no valid value for 'alpha'. For help see ?ppclust")
  
  sample.size <- apply(data, 1, function(x) length(x[!is.na(x)]))
  
  if(any(sample.size < 2))
    stop("all sample sizes must be larger than 2")
  
 
  if (!(alpha > 0 & alpha < 1))
    stop("'alpha' must be a real number in the range (0,1)")

  anovarank <- function(data, colms, cols4, colv22, coltal2)
  {
    if (is.vector(data)) return(1)
    sigma4 <- data[, cols4]
    a <- nrow(data)
    Ri. <- data[, colms]
    R.. <- mean(Ri.)
    MSTr <- 1 / (a - 1) * sum((Ri. - R..) ** 2)
    S2ri <- data[, colv22]
    MSEr <- 1 / a * sum(S2ri)
    Fr <- MSTr / MSEr
    est <- sqrt(a) * (Fr - 1)
    v22 <- 1 / a * sum(S2ri)
    tal2p <- data[, coltal2]
    tal2 <- 1 / a * sum(tal2p)
    var <- tal2 / v22 ** 2
    pv <- 1 - pnorm(est / sqrt(var))
    return(pv)
  }

  sortmed <- function(data, colms2)
  {
    data[order(data[, colms2]), ]
  }

  sigma4est <- function(x)
  {
    n <- length(x)
    if (n > 1) {
      s4hat <- (1 / n * sum((x - mean(x)) ** 2)) ** 2
      s4jack <- n * s4hat
      m <- matrix(c(rep(x, n)), byrow = T, ncol = n)
      diag(m) <- NA
      s4hatni <-
        apply(m, 1, function(x) {
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

  partition <- function(data, cts, cms)
  {
    if (nrow(data) < 2) {
      return(list(pt = 1, cts = cts))
    } else{
      datax <- data[, cms]
      count <- list(cts, cts + 1)
      n.int <- 0
      while (!all(count[[1]] %in% count[[2]]) & n.int < 13) {
        environment(cts) <- .GlobalEnv
        clusterExport(cl, varlist = c("cts"), envir = environment())
        mst.p <- parApply(cl, matrix(datax, ncol = 1), 1,
                                    function(x) (as.numeric(cts) - x) ** 2)
        mst.p <- t(mst.p)

        pt <- apply(mst.p, 1, function(x) which.min(x))
        cts <- tapply(datax, pt, function(x) median(x))
        for (i in 1:length(cts)) cts[i] <- data[which.min(abs(datax - cts[i])), cms]
        cts <- matrix(cts, ncol = 1)
        count[[1]] <- count[[2]]
        count[[2]] <- cts
        n.int <- n.int + 1
      }
      return(list(pt = as.numeric(pt), cts = cts))
    }
  }

  groupopt <- function(data, colgr, alpha)
  {
    n <- length(unique(data[, colgr]))
    m <- matrix(0, nrow = n, ncol = n)
    indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
    if (nrow(indices) == 1){
      indices <- matrix(indices[order(indices[, 2], indices[, 1]), ], ncol = 2)
    } else{
      indices <- indices[order(indices[, 2], indices[, 1]), ]
    }
    pvs <- apply(indices, 1, function(x)
      anovarank(data[data[, colgr] == x[1] |
                       data[, colgr] == x[2], ], cms, cs4, cv22, ctal2))
    m[lower.tri(m)] <- pvs
    mpv <- max(m, na.rm = T)
    while (mpv >= alpha) {
      idxs <- c(which(m == max(m, na.rm = T), arr.ind = T))
      x <- data[, colgr]
      x[x == max(idxs)] <- n + 1
      x[x == min(idxs)] <- n + 1
      x <- as.numeric(factor(x))
      data[, colgr] <- x
      m <- m[-idxs,-idxs]

      n <- n - 1
      indices <- matrix()
      indices <- cbind(rep(n, (n - 1)), 1:(n - 1))
      pvs <- apply(indices, 1, function(x)
        anovarank(data[data[, colgr] == x[1] |
                         data[, colgr] == x[2], ], cms, cs4, cv22, ctal2))
      m <- try(rbind(m, pvs), silent = T)
      m <- cbind(m, numeric(nrow(m)))
      mpv <- max(m, na.rm = T)
    }

    n <- length(unique(data[, colgr]))
    m <- matrix(NA, nrow = n, ncol = n)
    indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
    pvs <- apply(indices, 1, function(x)
      anovarank(data[data[,colgr]==x[1]|data[,colgr]==x[2],], cms, cs4, cv22, ctal2))
    m[lower.tri(m)] <- pvs
    indices <- cbind(1:n, 1:n)
    pvs <- apply(indices, 1, function(x)
      anovarank(data[data[,colgr]==x[1]|data[,colgr]==x[2],], cms, cs4, cv22, ctal2))
    diag(m) <- pvs
    m <- as.data.frame(matrix(format.pval(c(m), digits = 6),
                              nrow = n, ncol = n))
    colnames(m) <- 1:n
    list(cluster = data[, cgr], pv.matrix = m)
  }

  fib <- function(fibseq, x)
  {
    c(fibseq, fibseq[x - 1] + fibseq[x - 2])
  }

  cl <- makeCluster(n.cores)
  rankData <- matrix(rank(c(as.matrix(data)), na.last = 'keep'),
                     nrow = nrow(data), ncol = ncol(data))
  ncole <- ncol(rankData)
  nrowe <- nrow(rankData)
  environment(computeSigma4) <- .GlobalEnv
  clusterExport(cl, varlist = c("computeSigma4"), envir = environment())
  
  X <- parApply(cl, rankData, 1, function(x) {
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
  
  ppData <- data.table(id = 1:nrowe,
                       mean = means,
                       sigma4 = sigma4,
                       v22p = v22p,
                       tal2p = tal2p,
                       gr = 99999)
  # cid <- 1
  # cms <- 2
  # cs4 <- 3
  # cv22 <- 4
  # ctal2 <- 5
  # cgr <- 6
  fstart <- 1
  fend <- nrowe
  g <- 1
  
  pValue <- anovaRank(ppData)
  
  if (pValue > alpha) {
    stopCluster(cl)
    ppData[, "gr"] <- g
    return(list(cluster = ppData[, gr],
                P.values = pValue))
  } else{
    
    ppData <- sortMean(ppData)
    mid <- ppData[c(1, fend), mean]
    new.mid <- c()
    k <- 3
    fibseq <- c(1, 1)
    for (i in 1:15) fibseq <- fib(fibseq, i)
    
    cat('clustering ... \n')
    pb <- txtProgressBar(min = 0, max = fend, style = 3)
    while (fstart <= fend) {
      mid <- matrix(c(mid, new.mid), ncol = 1)
      groups <- partition(ppData[fstart:fend, ], mid, cms)
      ptt <- groups$pt
      mid <- groups$cts
      mg <- 0
      pValue <- numeric()
      
      for (i in 1:length(unique(ptt))) {
        pValue[i] <- anovaRank(ppData[c(fstart - 1 + which(ptt == i)), ],
                               cms, cs4, cv22, ctal2)
        if (pValue[i] > alpha**(1/(2*pi)) | is.nan(pValue[i])) {
          ppData[c(fstart - 1 + which(ptt == i)), gr] <- g
          g <- g + 1
          ppData <- ppData[order(ppData[, gr]), ]
          mg <- 1
        }
      }
      fstart <- sum(ppData[, gr] != 99999) + 1
      if (mg != 0) {
        mid <- c()
        k <- 3
      } else{
        mid <- c()
        k <- k + 1
      }
      if (fstart < fend)
        new.mid <- ppData[fstart:(fstart - 1 + fibseq[k]), ms]
      setTxtProgressBar(pb, fstart)
    }
    setTxtProgressBar(pb, fend)
    close(pb)
  }
  
  stopCluster(cl)
  
  resul <- groupOpt(ppData, gr, alpha)
  ppData[, gr] <- resul$cluster
  
  ppData <- ppData[order(ppData[, gr]), ]
  groupclass <- ppData[, c(cid, cms, gr)]
  groupclass0 <- groupclass[groupclass[, 3] == 0, ]
  groupclass <- groupclass[groupclass[, 3] != 0, ]
  groupc <- tapply(groupclass[, 2], groupclass[, 3], function(x) mean(x))
  groupc <- as.numeric(names(sort(groupc)))
  leng <- length(groupc)
  for(i in 1:leng){
    groupclass[groupclass[, 3] == groupc[i], ] <- leng + i
  }
  
  ppData[, gr] <- c(groupclass0[, 3], as.numeric(as.factor(groupclass[, 3])))
  ppData <- ppData[order(ppData[, cid]), ]
  
  return(list(cluster = ppData[, gr],
              P.values = resul$pv.matrix))
}

# Internal functions ------------------------------------------------------------------------

anovaRank <- function(x)
{
  sigma4 <- x[, "sigma4"]
  a <- nrow(x)
  Ri. <- x[, "mean"]
  R.. <- mean(Ri.)
  MSTr <- 1 / (a - 1) * sum((Ri. - R..) ** 2)
  S2ri <- x[, "v22p"]
  MSEr <- 1 / a * sum(S2ri)
  Fr <- MSTr / MSEr
  est <- sqrt(a) * (Fr - 1)
  v22 <- 1 / a * sum(S2ri)
  tal2p <- x[, "tal2p"]
  tal2 <- 1 / a * sum(tal2p)
  var <- tal2 / v22 ** 2
  pv <- 1 - pnorm(est / sqrt(var))
  return(pv)
}

sortMean <- function(x)
{
  x[order(x[, "mean"]), ]
}

computeSigma4 <- function(x)
{
  n <- length(x)
  if (n > 1) {
    s4hat <- (1 / n * sum((x - mean(x)) ** 2)) ** 2
    s4jack <- n * s4hat
    m <- matrix(c(rep(x, n)), byrow = T, ncol = n)
    diag(m) <- NA
    s4hatni <-
      apply(m, 1, function(x) {
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

partition <- function(data, cts, cms)
{
  if (nrow(x) < 2) {
    return(list(pt = 1, cts = cts))
  } else{
    datax <- x[, mean]
    count <- list(cts, cts + 1)
    n.int <- 0
    while (!all(count[[1]] %in% count[[2]]) & n.int < 13) {
      environment(cts) <- .GlobalEnv
      clusterExport(cl, varlist = c("cts"), envir = environment())
      mst.p <- parApply(cl, matrix(datax, ncol = 1), 1,
                        function(x) (as.numeric(cts) - x) ** 2)
      mst.p <- t(mst.p)
      
      pt <- apply(mst.p, 1, function(x) which.min(x))
      cts <- tapply(datax, pt, function(x) median(x))
      for (i in 1:length(cts)) cts[i] <- x[which.min(abs(datax - cts[i])), mean]
      cts <- matrix(cts, ncol = 1)
      count[[1]] <- count[[2]]
      count[[2]] <- cts
      n.int <- n.int + 1
    }
    return(list(pt = as.numeric(pt), cts = cts))
  }
}

groupOpt <- function(data, colgr, alpha)
{
  n <- length(unique(data[, gr]))
  m <- matrix(0, nrow = n, ncol = n)
  indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
  if (nrow(indices) == 1){
    indices <- matrix(indices[order(indices[, 2], indices[, 1]), ], ncol = 2)
  } else{
    indices <- indices[order(indices[, 2], indices[, 1]), ]
  }
  pvs <- apply(indices, 1, function(x) anovaRank(data[data[, gr] == x[1] | data[, gr] == x[2], ]))
  m[lower.tri(m)] <- pvs
  mpv <- max(m, na.rm = T)
  while (mpv >= alpha) {
    idxs <- c(which(m == max(m, na.rm = T), arr.ind = T))
    x <- data[, gr]
    x[x == max(idxs)] <- n + 1
    x[x == min(idxs)] <- n + 1
    x <- as.numeric(factor(x))
    data[, gr] <- x
    m <- m[-idxs,-idxs]
    
    n <- n - 1
    indices <- matrix()
    indices <- cbind(rep(n, (n - 1)), 1:(n - 1))
    pvs <- apply(indices, 1, function(x) anovaRank(data[data[, gr] == x[1] | data[, gr] == x[2], ]))
    m <- try(rbind(m, pvs), silent = T)
    m <- cbind(m, numeric(nrow(m)))
    mpv <- max(m, na.rm = T)
  }
  
  n <- length(unique(data[, gr]))
  m <- matrix(NA, nrow = n, ncol = n)
  indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
  pvs <- apply(indices, 1, function(x) anovaRank(data[data[, gr] == x[1] | data[, gr] == x[2], ]))
  m[lower.tri(m)] <- pvs
  indices <- cbind(1:n, 1:n)
  pvs <- apply(indices, 1, function(x) anovaRank(data[data[, gr] == x[1] | data[, gr] == x[2], ]))
  diag(m) <- pvs
  m <- as.data.frame(matrix(format.pval(c(m), digits = 6),
                            nrow = n, ncol = n))
  colnames(m) <- 1:n
  list(cluster = data[, gr], pv.matrix = m)
}

fib <- function(fibseq, x)
{
  c(fibseq, fibseq[x - 1] + fibseq[x - 2])
}

