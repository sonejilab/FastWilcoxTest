// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
#include <math.h>
using namespace Rcpp;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;


//' @name ZScore
//' @aliases ZScore,FastWilcoxTest-method
//' @rdname ZScore-methods
//' @docType methods
//' @description A specific z. score method that converts the data to 10 +-1 instead of 0+-1
//' in order to keep the not expressed clearly separate from the real data.
//' @param data the sparse Matrix
//' @param display_progress show a progress bar (TRUE)
//' @title Calculate z score for a sparse matrix
//' @export
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ZScore (Eigen::SparseMatrix<double> data, bool display_progress=true){
	/* Progress p(data.outerSize(), display_progress); */
	data = data.transpose();
	for (int k=0; k < data.outerSize(); ++k){
		/*p.increment();*/
		double sum = 0.0;
		int c = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( it.value() > 0) {
				c++;
				sum += it.value();
			}
		}
		double mean = sum / c;
		sum = 0.0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( it.value() > 0) {
				double entry = (it.value() - mean);
				sum += entry * entry;
				it.valueRef() = entry;
			}
		}
		double sd = sqrt(sum/c);
		/*Rcout << k << " mean " << mean << " and sd " << sd << "with count "<< c<< std::endl;*/
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( it.value() > 0) {
				double entry = (it.value() / sd) +10.0 ;
				it.valueRef() = entry;
			}
		}
	}
	data = data.transpose();
    return (data);
}
