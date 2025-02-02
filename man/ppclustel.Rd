% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppclustel.R
\name{ppclustel}
\alias{ppclustel}
\title{Clustering Algorithm to HDLLSS data.}
\usage{
ppclustel(dataset, id, rep, alpha, ...)
}
\arguments{
\item{dataset}{A numeric matrix or data frame with all numeric columns. If a matrix or data frame, rows correspond to variables (d) and columns correspond to observations (n).}

\item{id}{A integer number specifying the column of dataset with the variable id.}

\item{rep}{A integer number specifying the column of dataset indicating each variable replication number.}

\item{alpha}{A real number in the range (0, 1) indicanting the threshold parameter to be compared with p-values in the clustering procedure.}

\item{...}{not used.}
}
\value{
Results

A vector of integers indicating the cluster to which each variable is allocated and a cluster cross table of P-values.
}
\description{
Clustering Algorithm to HDLLSS data.
}
