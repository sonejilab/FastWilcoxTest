\name{logFC}
\alias{logFC}
\title{Fast log fold change}
\usage{logFC(A, B)}
\arguments{
  \item{A}{Numeric vector representing the population of interest.}
  \item{B}{Numeric vector representing the background population.}
}
\value{Numeric natural-log fold change.}
\description{Historical FastWilcoxTest convention: \code{log(mean(B) / mean(A))}.}
