#' Clustering Algorithm to HDLLSS data.
#'
#' @param dataset A matrix of data.
#' @param id A number.
#' @param rep A number.
#' @param alpha Significance (alpha) level.
#' @return A vector of integers indicating the cluster to which each point is allocated.
#' @return A cross table of P-values.

ppclustel <- function(dataset, id, rep, alpha) {

  anovalong <- function(resp, b)
  {
    if (is.vector(resp)) return(1)
    a <- nrow(resp)
    cts <- 5
    cte <- 4 + b
    xijp <- resp[, cts:cte]
    xpjp <- c(apply(xijp, 2, function(x) mean(x, na.rm = T)))
    mst <- c(apply(resp[, cts:cte], 1, function(x) (x - xpjp) ** 2))
    mst <- sum(mst) / ((a - 1) * b)
    mse <- sum(resp[, 3]) / a
    cov2 <- resp[, 4]
    cov2 <- sum(cov2) / a
    pv <- 1 - pnorm(sqrt(a * b) * (mst - mse) / sqrt(cov2))
    return(pv)
  }

  sortmed <- function(data)
  {
    data[order(data[, 2]),]
  }

  sortcenter <- function(data, tst)
  {
    nr <- nrow(data)
    nc <- ncol(data)
    if (tst == 1) {
      m2 <- data
    } else{
      ed <- tst-1
      m1 <- data[1:ed, ]
      m2 <- data[tst:nr, ]
    }
    nr2 <- nrow(m2)
    j1 <- numeric(nr2)
    c1 <- round(0.35 * nr2, 0)
    c2 <- round(0.65 * nr2, 0)
    j1[c(1:nr2) < c1] <- 1
    j1[c(1:nr2) > c2] <- 1
    m2 <- m2[order(j1, m2[, 2]), ]
    if (tst == 1)
      return(m2)
    else
      return(rbind(m1, m2))
  }

  indtest <- function(data, b, st, alpha, a, colgr, colts, g)
  {
    data[, colts] <- 0
    j1 <- data[, colts]
    j1[data[, colgr] == g] <- 1
    data[, colts] <- j1
    count <- 0

    if (st > a) return(list(st = st, data = data))

    for (i in st:a) {
      data[i, colts] <- 1
      data2 <- data[data[, colts] != 0, ]
      pvalue <- do.call('anovalong', list(data2, b))
      if (pvalue > alpha) data[i, colgr] <- g
      else data[i, colts] <- 0
    }
    data <- data[order(data[, colgr]), ]
    count <- sum(data[, colgr] <= g)
    st <- count + 1
    list(st = st, data = data)
  }

  t_start <- 3
  t_end <- ncol(dataset)
  names(dataset)[c(id, rep)] <- c('gene', 'array')
  dataset <- as.data.frame(dataset)
  nigene <-
    data.frame(gene = unique(dataset[, 1]), count = as.numeric(table(dataset[, 1])))
  names(dataset)[t_start:t_end] <- paste('time.', seq_along(t_start:t_end), sep = '')
  datasetf2 <- reshape(dataset, direction = "long", varying = t_start:t_end,
                       v.names = "exp", timevar = "time")
  datasetf2 <- datasetf2[, -5]
  b <- length(t_start:t_end)

  result1 <- aggregate(datasetf2[, 4], list(datasetf2[, id]),
                       function(x) median(x, na.rm = T))
  names(result1) <- c('gene','median')
  result2 <- aggregate(datasetf2[, 4], list(datasetf2[, id], datasetf2[, 3]),
                       function(x) mean(x, na.rm = T))
  result2 <- result2[order(result2[,1], result2[,2]),]
  names(result2) <- c('gene', 'time', 'rbijp')

  temp1 <- merge(datasetf2, result2, by = c('gene', 'time'))
  temp1 <- temp1[,c(1, 3, 2, 4:5)]
  temp1 <- temp1[order(temp1[,1], temp1[,3], temp1[,2]),]

  temp2 <- merge(temp1, nigene, by = 'gene')
  temp2$msepartial <- with(temp2, (exp - rbijp) ** 2 / (b * count * (count - 1)))

  result3 <- aggregate(temp2[, "msepartial"], list(temp2[, "gene"]),
                       function(x) sum(x, na.rm = T))
  names(result3) <- c('gene', 'msep')

  a <- length(unique(dataset[, 1]))
  cov2 <- c()
  for (i in 1:a) {
    cv <- cov(dataset[dataset[, 1] == i, t_start:t_end])
    cov2 <- c(cov2, sum(cv ** 2))
  }

  count <- nigene$count
  result4 <-
    data.frame(gene = unique(dataset[, 1]), cov2p = 2 * cov2 / (b * count * (count - 1)))
  result5 <-
    aggregate(dataset[, t_start:t_end], list(dataset[, id]),
              function(x) mean(x, na.rm = T))
  names(result5)[1] <- "gene"

  mdata <- merge(result1, result3, by = 'gene')
  mdata <- merge(mdata, result4, by = 'gene')
  mdata <- merge(mdata, result5, by = 'gene')

  g <- 1
  a <- nrow(mdata)
  ncole <- ncol(mdata)
  nrowe <- nrow(mdata)
  mdata <- do.call('sortmed', list(mdata))
  colid <- 1
  colgr <- ncole + 1
  colts <- ncole + 2
  mdata <- data.frame(mdata, gr = rep(99999, a), ts = rep(0, a))
  tstart <- 1
  tend <- a
  ktest <- 0
  mdata <- do.call('sortcenter', list(mdata, tstart))
  pv <- do.call('anovalong', list(mdata, b))
  pv2 <- pv
  if (pv > alpha) {
    mdata[, colgr] <- g
    return(list(cluster = mdata[, colgr],
                P.values = pv))
  } else{

    cat('clustering ... \n')
    L <- (tend - tstart + 1) / 2
    tend <- tstart + floor(L - 1)
    pb <- txtProgressBar(min = 0, max = a, style = 3)
    while (tstart <= a){
      if(tstart == a) {
        mdata[tstart, colgr] <- 0
        tstart <- tstart + 1
      } else{
        if (pv2 > alpha) g <- g + 1
        mdata <- do.call('sortcenter', list(mdata, tstart))
        mdata[, colts] <- 0
        j1 <- mdata[, colts]
        j1[which(tstart <= c(1:a) & c(1:a) <= tend)] <- 1
        mdata[, colts] <- j1

        n <- sum(mdata[, colts])
        if (n == 1){
          mdata[tstart, colgr] <- 0
          tstart <- tend + 1
          tend <- a
          ktest <- -1
        }

        exp <- mdata[mdata[, colts] != 0, ]
        if (ktest == 0) pv <- do.call('anovalong',list(exp, b))
        pv2 <- pv

        if (pv > alpha){
          j2 <- mdata[, colgr]
          j2[which(tstart <= c(1:a) & c(1:a) <= tend)] <- g
          mdata[, colgr] <- j2
          tstart <- tend + 1
          tend <- a
          setTxtProgressBar(pb, tstart)
          resul <- do.call('indtest', list(mdata, b, tstart, alpha, a, colgr, colts, g))
          tstart <- resul$st
          mdata <- resul$data
          setTxtProgressBar(pb, tstart)
        } else{
          L <- (tend - tstart + 1) * 0.9
          tend <- tstart + floor(L - 1)
        }
      }
    }
  }

  setTxtProgressBar(pb, a)
  close(pb)

  mdata <- mdata[order(mdata[, colgr]), ]
  groupclass <- mdata[, c(1, 2, colgr)]
  groupclass0 <- groupclass[groupclass$gr == 0, ]
  groupclass <- groupclass[groupclass$gr != 0, ]
  groupc <- tapply(groupclass$median, groupclass$gr, function(x) mean(x))
  groupc <- as.numeric(names(sort(groupc)))
  leng <- length(groupc)
  for(i in 1:leng){
    groupclass[groupclass$gr == groupc[i], ] <- leng + i
  }

  mdata[, colgr] <- c(groupclass0$gr, as.numeric(as.factor(groupclass$gr)))
  mdata <- mdata[order(mdata[, colid]), ]

  n <- length(unique(mdata[, colgr]))
  st <- min(unique(mdata[, colgr]))
  en <- max(unique(mdata[, colgr]))
  m <- matrix(NA, nrow = n, ncol = n)
  indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1))) + (st - 1)
  pvs <- apply(indices, 1, function(x)
    anovalong(mdata[mdata[, colgr] == x[1] |
                      mdata [, colgr] == x[2], ], b))
  m[lower.tri(m)] <- pvs
  indices <- cbind(st:en, st:en)
  pvs <- apply(indices, 1, function(x)
    anovalong(mdata[mdata[, colgr] == x[1] |
                      mdata [, colgr] == x[2], ], b))
  diag(m) <- pvs
  m <- as.data.frame(matrix(format.pval(c(m), digits = 6),
                            nrow = n, ncol = n))
  colnames(m) <- st:en
  row.names(m) <- colnames(m)

  return(list(cluster = mdata[order(mdata[, colid]), colgr],
              P.values = m))
}
