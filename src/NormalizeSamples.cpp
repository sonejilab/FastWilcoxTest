#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>
using namespace Rcpp;
#include <progress.hpp>
//#include <Rcpp::sugar::sample.h>



//' @name NormalizeSamples
//' @aliases NormalizeSamples,FastWilcoxTest-method
//' @rdname NormalizeSamples-methods
//' @docType methods
//' @description Normalize the NGS expression values by log re-scaling.
//' @param data the sparse Matrix (row = genes, col = cells)
//' @param scaleFactor the geometric mean of the data (use DESeq2 to get these)
//' @param display_progress show a progress bar (TRUE)
//' @title rescale a matrix using scaleFactor per column
//' @export
//[[Rcpp::export]]
Eigen::SparseMatrix<double>  NormalizeSamples (Eigen::SparseMatrix<double> X, std::vector<double> scaleFactor, bool display_progress=true ){

	Progress p(X.outerSize(), display_progress);
	
	for (int k=0; k < X.outerSize(); ++k){

		p.increment();
		
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it){
			double neu;
			neu = it.value() / scaleFactor[k];
			it.valueRef() = neu;
		}
	}
	return (X);
}
		