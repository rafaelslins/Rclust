\name{ppclustel_h}
\alias{ppclustel_h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Algorithm based on P-values.
}
\description{
The function \code{ppclustel_h()} perform clustering on a high dimensional data matrix with repeated measuriments on time.
}
\usage{
ppclustel_h(dataset, id, rep, alpha, n.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
numeric matrix of data, or an object that can be coerced to such a matrix.}
 \item{id}{
the dataset col number that contains the genes id.}
  \item{rep}{
the dataset col number that identify the repeated measuriments.}
  \item{alpha}{
significance level to clustering.}
}
\details{
The \code{alpha} value belongs the open interval (0,1).

Put the test statistic...
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
von Borries G. (2008) Partition clustering of High Dimensional Low Sample Size data based on P-Values. PhD Thesis.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# nf <- 1500
# b <- 10
# ni <- 3
# s1 <- matrix(rep(1, b * b), ncol = b)
# s2 <- exp(-1 / b * abs(matrix(rep(1:b, b), byrow = T, ncol = b) -
#                          matrix(rep(1:b, b), byrow = F, ncol = b)))
# Sig <- s1 * s2
# m <- b
# j <- 1:m
# tij <- j/m
# mu1 = c()
# for(j in tij)
#   mu1 <- c(mu1, min((2 - 5 * j) / 2, ((5 * j - 2) / 3) ^ 2 + sin((5 * pi *j) / 2)))
# mu2 <- -mu1
# mu3 <- cos(2 * pi * tij)
# mu4 <- -mu3
# g1 <- mvrnorm(n = 0.25 * nf * ni, mu1, Sig)
# g2 <- mvrnorm(n = 0.25 * nf * ni, mu2, Sig)
# g3 <- mvrnorm(n = 0.25 * nf * ni, mu3, Sig)
# g4 <- mvrnorm(n = 0.25 * nf * ni, mu4, Sig)
# fac <- rep(1:nf, each = nobs)
# rep <- rep(1:ni, nf)
# simullongdata <- cbind(fac, rep, rbind(g1, g2, g3, g4))
#
# cl <- ppclustel_h(simullongdata, id = 1, rep = 2, alpha = 1e-5, n.cores = 1)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}
