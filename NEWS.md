# FastWilcoxTest 0.3.0

- Modernized the sparse Wilcoxon core for current R/Matrix stacks.
- Replaced RcppEigen sparse mapping with direct `dgCMatrix` CSC access via Rcpp.
- Preserved the historical `StatTest()` API and fold-change convention.
- Removed legacy compiled artifacts and mixed C++ standard settings.
- Reduced the core dependency set.
- Added regression tests against `stats::wilcox.test(exact = FALSE)`.
