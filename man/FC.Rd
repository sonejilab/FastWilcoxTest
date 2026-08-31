\name{FC}
\alias{FC}
\title{Fast fold change}
\usage{FC(A, B)}
\arguments{
  \item{A}{Numeric vector representing the population of interest.}
  \item{B}{Numeric vector representing the background population.}
}
\value{Numeric fold change, \code{mean(B) / mean(A)}.}
\description{Historical FastWilcoxTest fold-change convention: mean(B) divided by mean(A).}
