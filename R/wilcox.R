#' Fast fold change
#'
#' Historical FastWilcoxTest fold-change convention: mean(B) / mean(A).
#'
#' @param A Numeric vector representing the population of interest.
#' @param B Numeric vector representing the background population.
#' @return Numeric fold change, `mean(B) / mean(A)`.
#' @export
FC <- function(A, B) {
  .Call("_FastWilcoxTest_fc", as.numeric(A), as.numeric(B), PACKAGE = "FastWilcoxTest")
}

#' Fast log fold change
#'
#' Historical FastWilcoxTest log-fold-change convention:
#' `log(mean(B) / mean(A))` using the natural logarithm.
#'
#' @inheritParams FC
#' @return Numeric natural-log fold change.
#' @export
logFC <- function(A, B) {
  .Call("_FastWilcoxTest_log_fc", as.numeric(A), as.numeric(B), PACKAGE = "FastWilcoxTest")
}

#' Fast Wilcoxon rank-sum test for two vectors
#'
#' @param x Numeric vector for the population of interest.
#' @param y Numeric vector for the background population.
#' @param type Integer test result type. `0` = greater, `1` = less,
#'   `2` = two-sided, `3` = U statistic, `4` = -log10(greater),
#'   `5` = log10(less), `6` = -log10(two-sided), `7` = signed
#'   -log10(two-sided).
#' @return Numeric vector containing rank sum and requested statistic/p-value.
#' @export
cppWilcoxTest <- function(x, y, type = 2L) {
  .Call(
    "_FastWilcoxTest_cpp_wilcox_test",
    as.numeric(x),
    as.numeric(y),
    as.integer(type),
    PACKAGE = "FastWilcoxTest"
  )
}

#' Fast Wilcoxon marker tests on a sparse matrix
#'
#' Performs Wilcoxon rank-sum tests on selected columns of a sparse matrix.
#' Rows represent observations/cells and columns represent features/genes.
#' Expression-frequency and fold-change filters are evaluated before rank tests
#' are run.
#'
#' @param X Sparse numeric matrix with observations in rows and features in
#'   columns. Objects inheriting from `Matrix::sparseMatrix` are coerced to
#'   `dgCMatrix` without densifying.
#' @param interest 1-based row indices for the population of interest.
#' @param background 1-based row indices for the background population.
#' @param logFCcut Minimum absolute natural-log fold change required for a
#'   feature to be tested. Historical default is `1`.
#' @param minPct Minimum detected fraction in either population. A value is
#'   detected when it is greater than zero.
#' @param onlyPos If `TRUE`, retain only features satisfying the historical
#'   positive direction (`mean(background) / mean(interest)` above the cutoff).
#'
#' @return Numeric matrix with columns `colID`, `logFC`, `fracExprIN`,
#'   `fracExprOUT`, `rank.sum`, `p.value`, and `FC`. `colID` is 1-based.
#'
#' @details
#' Version 0.3 reads the compressed sparse column representation directly and
#' no longer maps Matrix objects through RcppEigen. This keeps the calculation
#' sparse and removes a fragile dependency boundary while preserving the public
#' API of the original implementation.
#'
#' For compatibility, fold change keeps the original FastWilcoxTest convention:
#' `FC = mean(background) / mean(interest)` and `logFC = log(FC)`.
#'
#' @examples
#' X <- Matrix::Matrix(
#'   c(10, 9, 8, 0, 0, 0,
#'     0, 0, 1, 8, 9, 10),
#'   nrow = 6,
#'   sparse = TRUE
#' )
#' StatTest(X, 1:3, 4:6, logFCcut = 0, minPct = 0)
#'
#' @export
StatTest <- function(
    X,
    interest,
    background,
    logFCcut = 1.0,
    minPct = 0.1,
    onlyPos = FALSE) {

  if (!inherits(X, "sparseMatrix")) {
    stop("X must be a sparse Matrix object", call. = FALSE)
  }

  if (!inherits(X, "dgCMatrix")) {
    X <- methods::as(X, "dgCMatrix")
  }

  interest <- as.integer(interest)
  background <- as.integer(background)

  if (!length(interest)) {
    stop("No values in interest group", call. = FALSE)
  }
  if (!length(background)) {
    stop("No values in background group", call. = FALSE)
  }
  if (anyNA(interest) || anyNA(background)) {
    stop("interest and background must not contain NA", call. = FALSE)
  }
  if (length(intersect(interest, background))) {
    stop("interest and background must not overlap", call. = FALSE)
  }
  if (length(logFCcut) != 1L || !is.finite(logFCcut) || logFCcut < 0) {
    stop("logFCcut must be one finite non-negative number", call. = FALSE)
  }
  if (length(minPct) != 1L || !is.finite(minPct) || minPct < 0 || minPct > 1) {
    stop("minPct must be between 0 and 1", call. = FALSE)
  }

  .Call(
    "_FastWilcoxTest_stat_test",
    X,
    interest,
    background,
    as.numeric(logFCcut),
    as.numeric(minPct),
    as.logical(onlyPos),
    PACKAGE = "FastWilcoxTest"
  )
}
