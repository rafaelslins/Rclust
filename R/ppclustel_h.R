#' Clustering Algorithm to HDLLSS data.
#'
#' @param dataset A numeric matrix or data frame with all numeric columns. If a matrix or data frame, rows correspond to variables (d) and columns correspond to observations (n).
#' @param id A integer number specifying the column of dataset with the variable id.
#' @param rep A integer number specifying the column of dataset indicating each variable replication number.
#' @param alpha A real number in the range (0, 1) indicanting the threshold parameter to be compared with p-values in the clustering procedure.
#' @param n.cores A number processor cores (see detectCores).
#' @param ... not used.
#' 
#' @return Results
#' 
#' A vector of integers indicating the cluster to which each variable is allocated and a cluster cross table of P-values.
#' 
#' @importFrom stats aggregate cov median pnorm reshape var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel clusterExport makeCluster parApply parLapply stopCluster
#' @export

ppclustel_h <- function(dataset, id, rep, alpha, n.cores = 1, ...){

  if (any(apply(dataset, 2, is.numeric) == FALSE))
    stop('dataset contains non-numeric values')
  
  if (!(alpha > 0 & alpha < 1))
    stop("'alpha' must be a real number in the range (0,1)")
  
  anovalong <- function(data, t1, tb)
  {
    if (is.vector(data)) return(1)
    a <- nrow(data)
    b <- length(t1:tb)
    xpjp <- as.numeric(apply(data[, t1:tb], 2, function(x) mean(x, na.rm = T)))
    msp <- c(apply(data[, t1:tb],1, function(x) (x - xpjp)**2))
    msphi <- sum(msp) / ((a - 1) * b)
    mse <- sum(data[, 3]) / a
    cov2 <- sum(data[, 4]) / a
    pv <- 1 - pnorm(sqrt(a * b) * (msphi - mse) / sqrt(cov2))
    return(pv)
  }

  groupopt <- function(data, colgr, t1, tb, alpha)
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
      anovalong(data[data[, colgr] == x[1] | data[, colgr] == x[2], ], t1, tb))

    m[lower.tri(m)] <- pvs
    mpv <- max(m, na.rm = T)
    while (mpv >= alpha) {
      idxs <- c(which(m == max(m, na.rm = T), arr.ind = T))[1:2]
      x <- data[, colgr]
      x[x == max(idxs)] <- n + 1
      x[x == min(idxs)] <- n + 1
      x <- as.numeric(factor(x))
      data[, colgr] <- x
      m <- m[-idxs, -idxs]

      n <- n - 1
      indices <- matrix()
      indices <- cbind(rep(n, (n - 1)), 1:(n - 1))
      pvs <- apply(indices, 1, function(x)
        anovalong(data[data[, colgr] == x[1] |
                         data[, colgr] == x[2], ], t1, tb))
      m <- try(rbind(m, pvs), silent = T)
      m <- cbind(m, numeric(nrow(m)))
      mpv <- max(m, na.rm = T)
    }

    n <- length(unique(data[, colgr]))
    m <- matrix(NA, nrow = n, ncol = n)
    indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1)))
    pvs <- apply(indices, 1, function(x)
      anovalong(data[data[, colgr] == x[1] | data[, colgr] == x[2], ], t1, tb))
    m[lower.tri(m)] <- pvs
    indices <- cbind(1:n, 1:n)
    pvs <- apply(indices, 1, function(x)
      anovalong(data[data[, colgr] == x[1] | data[, colgr] == x[2], ], t1, tb))
    diag(m) <- pvs
    m <- as.data.frame(matrix(format.pval(c(m), digits = 6),
                              nrow = n, ncol = n))
    colnames(m) <- 1:n
    return(list(cluster = data[, cgr], pv.matrix = m))
  }


  partition <- function(data, cts, t1, tb)
  {
    if (nrow(data) < 2) {
      return(list(pt = 1, cts = cts))
    } else{
      datax <- as.matrix(data[, t1:tb])
      count <- list(cts, cts + 1)
      n.int <- 0
      while (!all(count[[1]] %in% count[[2]]) & n.int < 13) {
        mst.p <- apply(cts, 1, function(x) {
          dif2 <- t(t(datax) - as.numeric(x)) ** 2
          dif2 %*% matrix(rep(1, ncol(dif2)))
        })
        pt <- apply(mst.p, 1, function(x) which.min(x))
        cts <- aggregate(datax, list(pt), function(x) mean(x))[, -1]
        for (i in 1:nrow(cts)) {
          a <- t(t(datax) - as.numeric(cts[i, ])) ** 2
          a <- which.min(a %*% matrix(rep(1, ncol(a))))
          cts[i, ] <- datax[a, ]
        }
        count[[1]] <- count[[2]]
        count[[2]] <- cts
        n.int <- n.int + 1
      }
      return(list(pt = as.numeric(pt), cts = cts))
    }
  }

  fib <- function(fibseq, x)
  {
    c(fibseq, fibseq[x-1] + fibseq[x-2])
  }

  t_start <- 3
  t_end <- ncol(dataset)

  ds1 <- as.data.frame(dataset[, c(id, rep, t_start:t_end)])
  colnames(ds1)[1:2] <- c('factor', 'rep')
  fac <- c(factor(ds1[, 'factor']))
  ds1[, 'factor'] <- fac

  a <- length(unique(fac))
  b <- length(t_start:t_end)
  t1 <- 3
  tb <- ncol(ds1)

  count <- as.numeric(table(fac))

  ds1 <- split(ds1[, t1:tb], fac)

  cl <- makeCluster(n.cores)
  cov2 <- parLapply(cl, ds1, function(x) sum(cov(x) ** 2))
  cov2 <- as.numeric(unlist(cov2))
  cov2p <- 2 * cov2 / (b * count * (count - 1))

  result1 <- parLapply(cl, ds1,
                                 function(x) apply(x, 2,
                                                   function(y) mean(y, na.rm = T)))
  result1 <- matrix(unlist(result1), byrow = T, ncol = b, nrow = a)
  means <- parApply(cl, result1, 1, function(x) mean(x,na.rm = T))

  result2 <- parLapply(cl, ds1,
                                 function(x) apply(x, 2,
                                                   function(y) sum((y-mean(y))**2, na.rm = T)))
  result2 <- matrix(unlist(result2), byrow = T, ncol = b, nrow = a)
  msep <- parApply(cl, result2, 1, function(x) sum(x,na.rm = T))
  msep <- msep / (b * count * (count - 1))

  ppmatrix <- data.frame(unique(fac), means, msep, cov2p, rep(99999, a), result1)
  colnames(ppmatrix) <- c('factor', 'means', 'msep', 'cov2p', 'gr',
                          paste('t', 1:b, sep = ''))
  ncole <- ncol(ppmatrix)
  cid <- 1
  cms <- 2
  cgr <- 5
  fend <- a
  g <- 1
  t1 <- 6
  tb <- ncole

  options(digits = 22)
  pvalue <- do.call('anovalong',list(ppmatrix, t1, tb))
print(pvalue)
  if (pvalue > alpha) {
    stopCluster(cl)
    ppmatrix[, cgr] <- g
    return(list(cluster = ppmatrix[, cgr],
                P.values = pvalue))
  } else{

    fstart <- 1
    ppmatrix <- ppmatrix[order(ppmatrix[, cms]), ]
    mid <- ppmatrix[c(1, a), t1:ncole]
    new.mid <- numeric()
    k <- 3
    fibseq <- c(1, 1)
    for (i in 1:25) fibseq <- fib(fibseq, i)

    cat('clustering ...\n')
    pb <- txtProgressBar(min = 0, max = fend, style = 3)
    while (fstart <= fend){
      mid <- rbind(mid, new.mid)
      groups <- do.call('partition', list(ppmatrix[fstart:fend, ], mid, t1, tb))
      ptt <- groups$pt
      mid <- groups$cts
      mg <- 0
      pvalue <- c()
      for (i in 1:length(unique(ptt))) {
        pvalue[i] <- do.call('anovalong',
                             list(ppmatrix[c(fstart - 1 + which(ptt == i)), ], t1, tb))
        if (pvalue[i] > alpha**(1/(2*pi)) | is.nan(pvalue[i])) {
          ppmatrix[c(fstart - 1 + which(ptt == i)), cgr] <- g
          g <- g + 1
          mg <- 1
          ppmatrix <- ppmatrix[order(ppmatrix[, cgr]), ]
        }
      }
      fstart <- sum(ppmatrix[, cgr] != 99999) + 1
      if (mg != 0) {
        mid <- c()
        k <- 3
      } else{
        mid <- c()
        k <- k + 1
      }
      new.mid <- ppmatrix[fstart:(fstart - 1 + fibseq[k]), t1:tb]
      setTxtProgressBar(pb, fstart)
    }
    setTxtProgressBar(pb, fend)
    close(pb)

    stopCluster(cl)

    ppmatrix[, cgr] <- as.numeric(factor(ppmatrix[, cgr]))
    resul <- do.call('groupopt',list(ppmatrix, cgr, t1, tb, alpha))
    ppmatrix[, cgr] <- resul$cluster
    ppmatrix <- ppmatrix[order(ppmatrix[, cid]),]

    return(list(cluster = ppmatrix[, cgr],
                P.values = resul$pv.matrix))
  }
}
