\name{StatTest}
\alias{StatTest}
\title{Fast Wilcoxon marker tests on a sparse matrix}
\usage{
StatTest(X, interest, background, logFCcut = 1, minPct = 0.1,
  onlyPos = FALSE)
}
\arguments{
  \item{X}{Sparse numeric matrix with observations/cells in rows and features/genes in columns.}
  \item{interest}{1-based row indices for the population of interest.}
  \item{background}{1-based row indices for the background population.}
  \item{logFCcut}{Minimum absolute natural-log fold-change threshold.}
  \item{minPct}{Minimum detected fraction in either population.}
  \item{onlyPos}{If TRUE, retain only the historical positive fold-change direction.}
}
\value{Numeric matrix with columns \code{colID}, \code{logFC}, \code{fracExprIN}, \code{fracExprOUT}, \code{rank.sum}, \code{p.value}, and \code{FC}.}
\description{Performs sparse marker-style Wilcoxon rank-sum tests with filtering before testing.}
\details{Version 0.3 reads the compressed sparse column representation directly and no longer maps Matrix objects through RcppEigen. For compatibility, fold change remains \code{mean(background) / mean(interest)}.}
