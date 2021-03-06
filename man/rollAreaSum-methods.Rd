% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\docType{methods}
\name{rollAreaSum}
\alias{rollAreaSum}
\alias{rollAreaSum,FastWilcoxTest-method}
\title{rolling sum over sparse matrix}
\usage{
rollAreaSum(X, location, funcID, size)
}
\arguments{
\item{X}{the sparse matrix}

\item{location}{the location for every row in the matrix}

\item{funcID}{two functions : 1 == sum; 2 == mean}

\item{size}{the size of the rolling window}
}
\description{
calculate a rolling sum of the rows
}
\details{
The numbers start at the first row and end at the last row having a full sized widow
}
