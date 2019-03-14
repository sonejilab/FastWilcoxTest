# FastWilcoxTest

This package defines several functions that work on sparse Matrix R objects and try to be as fast as possible.

The main function is a StatTest which performs wilcox signed rank tests directly in c++, giving exactly the same results as the wilcox.test R function, but ~10x faster. The method has been copied from the BioQC package, but is exported from this package, both in R as ```FastWilcoxTest::StatTest()``` as well as in Rcpp as the FastWilcoxTest.h header file.

In addition this package exports a cppWilcoxTest working on two vectors instead of a sparse matrix and a vector, a FastCor function and a very special ZScore function. The ZScore function projects the values > 0 into 10+-1 space. This is done to preserve the sparse quality of the matrix.

## Install

You need the R devtools package.

```
devtools::install_github('stela2502/FastWilcoxTest')
```


