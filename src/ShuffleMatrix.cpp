// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <Rdefines.h>

using namespace Rcpp;

//' @title reshuffle data based on a sparse matrix assuming max double the amount of entries not being zero
//' @aliases ShuffleMatrix,FastWilcoxTest-method
//' @rdname ShuffleMatrix
//' @description replacing the synthetic1 function of RFclust.SGE package 
//' @param X the sparse matrix (tests are applied to columns!)
//' @return a matrix with x, j and i avalues to be put into a new sparse matrix
//' @export
// [[Rcpp::export]]
 Eigen::SparseMatrix<double> ShuffleMatrix (Eigen::SparseMatrix<double> X) {

	Eigen::SparseMatrix<double> mat(X.innerSize(), X.outerSize());
	int S=0;
	if ( X.nonZeros() < X.innerSize() * X.outerSize() ) {
		//allow at max twice the non zeros??
		S = X.nonZeros() * 2;
	}else {
		S = X.nonZeros();
	}
	Rcout << "reserving " <<  S << " entries" << std::endl;

	mat.reserve( S );

	//std::vector<double> A(X.innerSize());
	NumericVector A(X.innerSize());
	//NumericVector Resampled(1);
	int Resampled = 0;
	int prob =1;
	int i;
	
	for (int c_=0; c_ < X.outerSize(); ++c_){
		//std::fill(A.begin(), A.end(), 0.0);
		A.fill(0.0);

		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		//Rcout << "calculate randoms " <<  A.length() << std::endl;
		//NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
		// and now we need to fill that into the sparse matrix..
		for ( int i=0; i < A.length(); i++){
			Resampled =  std::floor(R::runif(0,1) * A.length()) ; // Rcpp::sugar::SampleReplace(A, 1, R::runif(0,A.length()), false);
			if ( A[Resampled] != 0.0) {
				//Rcout << "insert " <<  i << ","<< c_ <<" value " << A[Resampled] << std::endl;
				mat.insert( i, c_ ) = A[Resampled];
			}
		}
	}
	return ( mat );
}