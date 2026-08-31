# FastWilcoxTest

Fast Wilcoxon rank-sum testing for sparse single-cell matrices.

Version 0.3 is a 2026 modernization of the original package. The main goal is
still the same: perform marker-style Wilcoxon tests without expanding a sparse
single-cell matrix or iterating through genes in R.

## What changed in 0.3

- `StatTest()` keeps the historical public API.
- Sparse `dgCMatrix` data are read directly from their CSC `p`, `i`, and `x`
  slots through Rcpp.
- RcppEigen and RcppProgress are no longer required by the Wilcoxon core.
- Old generated Rcpp interface files and stale compiled `.o`/`.so` artifacts
  are gone.
- The build uses a single C++17 configuration on Unix and Windows.
- Runtime dependencies are reduced to Rcpp, Matrix, methods, and stats.
- The Wilcoxon core has focused regression tests against `stats::wilcox.test()`.

## Install

```r
devtools::install_github("sonejilab/FastWilcoxTest")
```

or from a local source checkout:

```r
devtools::install(".")
```

## Sparse marker test

`StatTest()` expects cells/observations in rows and genes/features in columns:

```r
X <- Matrix::Matrix(expression_matrix, sparse = TRUE)
res <- FastWilcoxTest::StatTest(
    X,
    interest = which(group == "cluster_0"),
    background = which(group != "cluster_0"),
    logFCcut = 0.25,
    minPct = 0.1
)
```

The package intentionally preserves the original fold-change direction:
`FC = mean(background) / mean(interest)`. Downstream wrappers that want Seurat's
usual positive-is-interest convention should invert `FC` and negate `logFC`.
