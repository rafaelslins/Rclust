\name{ppclust_h}
\alias{ppclust_h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Algorithm based on P-values.
}
\description{
The function \code{ppclust_h()} perform clustering on a high dimensional data matrix with repeated measuriments.
}
\usage{
ppclust_h(dataset, alpha, n.cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
numeric matrix of data, or an object that can be coerced to such a matrix.}
  \item{alpha}{
significance level to clustering.}
  \item{n.cores}{
number of CPU cores to use.}
}
\details{
The \code{alpha} (significance level) belongs the open interval (0, 1).

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

# nf <- 15000
# nobs <- 3
# g1 <- 0.25 * rnorm(0.20 * nf * nobs) - 0.5
# g2 <- 0.25 * rnorm(0.15 * nf * nobs) - 0.3
# g3 <- 0.25 * rnorm(0.30 * nf * nobs)
# g4 <- 0.25 * rnorm(0.20 * nf * nobs) + 0.5
# g5 <- 0.25 * rnorm(0.15 * nf * nobs) + 1
# simuldata <- matrix(c(g1, g2, g3, g4, g5), byrow = T, ncol = nobs, nrow = nf)

# cl <- ppclust_h(simuldata, alpha = 1e-5, n.cores = 1)

# table(cl$cluster)
#
# cl1$P.values

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}
