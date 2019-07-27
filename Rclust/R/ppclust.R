#' Clustering Algorithm to HDLSS data..
#'
#' @param dataset A data matrix.
#' @param alpha Significance (alpha) level.
#' @return A vector of integers indicating the cluster to which each point is allocated.
#' @return A cross table of P-values.

ppclust <- function(dataset, alpha) {

  if (any(apply(dataset, 2, is.numeric) == FALSE))
    stop('dataset contains non-numeric values')

  anovarank <- function(data, ncole, cols4)
  {
    if (is.vector(data)) return(1)
    sigma4 <- data[, cols4]
    data2 <- data[, 1:ncole]
    nfac <- nrow(data2)
    meang <- rowMeans(data2, na.rm = T)
    ng <- apply(data2, 1, function(x) length(x[!is.na(x)]))
    varg <- apply(data2, 1, function(x) var(x, na.rm = T))
    meant <- sum(meang) / nfac
    mst <- sum((meang - meant) ** 2) / (nfac - 1)
    mse2 <- sum(varg / ng) / nfac
    fr2 <- mst / mse2
    tau2 <- sum(sigma4 / (ng * (ng - 1))) * 2 / nfac
    asyvar <- tau2 / mse2 ** 2
    pv <- 1 - pnorm(sqrt(nfac) * (fr2 - 1) / sqrt(asyvar))
    return(pv)
  }

  sortmed <- function(data, ncole)
  {
    data2 <- data[, 1:ncole]
    data2[is.na(data2) == T] <- nrow(data2) + 1000
    medians <- apply(data2, 1, function(x) median(x, na.rm = T))
    data2 <- cbind(data, medians)
    data2[order(medians), ]
  }

  sortcenter <- function(data, ncole, colms, tst)
  {
    nr <- nrow(data[, 1:ncole])
    nc <- ncol(data[, 1:ncole])
    if (tst == 1) {
      m2 <- data
    } else{
      ed <- tst - 1
      m1 <- data[1:ed, ]
      m2 <- data[tst:nr, ]
    }
    nr2 <- nrow(m2)
    j1 <- numeric(nr2)
    c1 <- round(0.40 * nr2, 0)
    c2 <- round(0.60 * nr2, 0)
    j1[c(1:nr2) < c1] <- 1
    j1[c(1:nr2) > c2] <- 1
    m2 <- m2[order(j1, data[, colms]), ]
    if (tst == 1)
      return(m2)
    else
      return(rbind(m1, m2))
  }

  sigma4est <- function(data, ncole)
  {
    sig4hat <- function(x, n) {
      (sum((x - mean(x)) ** 2) / n) ** 2
    }

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
    sigma4 <- apply(data[, 1:ncole], 1, est)
    return(sigma4)
  }

  indtest <- function(data, st, nrowe, ncole, colgr, colts, cols4, g, pvalue, alpha)
  {
    data[, colts] <- 0
    j1 <- data[, colts]
    j1[data[, colgr] == g] <- 1
    data[, colts] <- j1
    count <- 0

    if (st > nrowe)
      return(list(st = st, data = data))

    for (i in st:nrowe) {
      data[i, colts] <- 1
      n <- sum(data[, colts])
      data2 <- data[data[, colts] != 0, ]
      pvalue <- do.call('anovarank', list(data2, ncole, cols4))
      if (pvalue > alpha)
        data[i, colgr] <- g
      else
        data[i, colts] <- 0
    }
    data <- data[order(data[, colgr]), ]
    count <- sum(data[, colgr] <= g)
    st <- count + 1
    list(st = st, data = data)
  }

  expression <- matrix(rank(c(as.matrix(dataset)), na.last = 'keep'),
                       nrow = nrow(dataset), ncol = ncol(dataset))
  ncole <- ncol(expression)
  nrowe <- nrow(expression)
  colid <- ncole + 1
  colgr <- ncole + 2
  cols4 <- ncole + 3
  colts <- ncole + 4
  colms <- ncole + 5
  tstart <- 1
  tend <- nrowe
  g <- 1
  testpp <- 0
  sigma4 <- do.call('sigma4est',list(expression, ncole))
  expression <- cbind(expression, id = 1:nrowe, gr = 99999, sigma4, colts = 0)
  expression <- as.data.frame(expression)
  expression <- do.call('sortmed',list(expression, ncole))
  expression <- do.call('sortcenter',list(expression, ncole, colms, tstart))
  pvalue <- do.call('anovarank',list(expression[tstart:tend,], ncole, cols4))

  if (pvalue > alpha) {
    expression[, colgr] <- 1
    return(list(cluster = expression[, colgr],
                P.values = pvalue))
  } else{

    cat('clustering ... \n')
    L <- (tend - tstart + 1) / 2
    tend <- tstart + floor(L - 1)
    pb <- txtProgressBar(min = 0, max = nrowe, style = 3)
    while (tstart <= nrowe) {
      if (tstart == nrowe) {
        expression[tstart, colgr] <- 0
        tstart <- tstart + 1
      } else{
        if (pvalue > alpha) g <- g + 1

        expression[, colts] <- 0
        j1 <- expression[, colts]
        j1[which(tstart <= c(1:nrowe) & c(1:nrowe) <= tend)] <- 1
        expression[, colts] <- j1

        n <- sum(expression[, colts])

        if (n == 1){
          expression[tstart, colgr] <- 0
          tstart <- tend + 1
          tend <- nrowe
          testpp <- -1
        }

        exp <- expression[expression[, colts] != 0,]
        if (testpp == 0) pvalue <- do.call('anovarank', list(exp, ncole, cols4))

        if (pvalue > alpha) {
          j2 <- expression[, colgr]
          j2[which(tstart <= c(1:nrowe) & c(1:nrowe) <= tend)] <- g
          expression[, colgr] <- j2

          tstart <- tend + 1
          tend <- nrowe
          setTxtProgressBar(pb, tstart)
          resul <- do.call('indtest', list(expression, tstart, nrowe, ncole,
                                           colgr, colts, cols4, g, pvalue, alpha))
          tstart <- resul$st
          expression <- resul$data
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

  expression <- expression[order(expression[, colgr]), ]
  groupclass <- expression[, c(colid, colms, colgr)]
  groupclass0 <- groupclass[groupclass$gr == 0, ]
  groupclass <- groupclass[groupclass$gr != 0, ]
  groupc <- tapply(groupclass$medians, groupclass$gr, function(x) mean(x))
  groupc <- as.numeric(names(sort(groupc)))
  leng <- length(groupc)
  for(i in 1:leng){
    groupclass[groupclass$gr == groupc[i], ] <- leng + i
  }

  expression[, colgr] <- c(groupclass0$gr, as.numeric(as.factor(groupclass$gr)))
  expression <- expression[order(expression[, colid]), ]

  n <- length(unique(expression[, colgr]))
  st <- min(unique(expression[, colgr]))
  en <- max(unique(expression[, colgr]))
  m <- matrix(NA, nrow = n, ncol = n)
  indices <- cbind(rep(2:n, 1:(n - 1)), sequence(1:(n - 1))) + (st - 1)
  pvs <- apply(indices, 1, function(x)
    anovarank(expression[expression[, colgr] == x[1] |
                           expression[, colgr] == x[2],], ncole, cols4))
  m[lower.tri(m)] <- pvs
  indices <- cbind(st:en, st:en)
  pvs <- apply(indices, 1, function(x)
    anovarank(expression[expression[, colgr] == x[1] |
                           expression[, colgr] == x[2],], ncole, cols4))
  diag(m) <- pvs
  m <- as.data.frame(matrix(format.pval(c(m), digits = 6), nrow = n, ncol = n))
  colnames(m) <- st:en
  row.names(m) <- colnames(m)
  return(list(cluster = expression[, colgr],
              P.values = m))
}
