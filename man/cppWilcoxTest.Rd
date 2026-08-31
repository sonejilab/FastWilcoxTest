\name{cppWilcoxTest}
\alias{cppWilcoxTest}
\title{Fast Wilcoxon rank-sum test for two vectors}
\usage{cppWilcoxTest(x, y, type = 2L)}
\arguments{
  \item{x}{Numeric vector for the population of interest.}
  \item{y}{Numeric vector for the background population.}
  \item{type}{Integer result type from 0 to 7. Type 2 is the usual two-sided p-value.}
}
\value{Named numeric vector containing \code{rank.sum} and \code{p.value}.}
\description{Runs the FastWilcoxTest Wilcoxon rank-sum implementation on two numeric vectors.}
