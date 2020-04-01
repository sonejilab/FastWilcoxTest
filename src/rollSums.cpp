// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
#include "FastWilcoxTest_RcppExports.h"


using namespace Rcpp;
#include <R.h>
#include <Rdefines.h>
#include <numeric>



//' @name rollSum
//' @aliases rollSum,FastWilcoxTest-method
//' @rdname rollSum-methods
//' @docType methods
//' @description calculate a rolling sum of the rows
//' @param X the sparse matrix
//' @param n the size of the rolling window
//' @title rolling sum over sparse matrix
//' @export
// [[Rcpp::export]]
NumericMatrix  rollSum (Eigen::SparseMatrix<double> X, int n){

	X= X.transpose();

	//Rcout << "calculating " <<  X.outerSize() << " tests (columns) using "<< CMP.size()<< " resp. " << X.innerSize() << " values" << std::endl;

	if ( X.innerSize() < n )
		::Rf_error("Sorry, the total columns (ncol(X)) in the data are less than the width of the rolling window (m) (%d, %d)", X.innerSize(), n );
	//NumericMatrix ret( nrow, ncol );
	NumericMatrix ret( X.outerSize(), X.innerSize() -n+1 );
	std::vector<double> A(X.innerSize());


	//Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;
	for (int c_=0; c_ < X.outerSize(); ++c_){
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		//Rcout << "summing " <<  X.innerSize() << " values up "<< n << " and col " << c_ << " values" << std::endl;
		//now iterate over all data
		for ( int i=n; i < X.innerSize()+1; i++){
			ret(c_, i-n) = std::accumulate(A.begin()+(i-n), A.begin()+i, 0.0);
			//Rcout << "sum from " <<  (i-n) <<" to "<< i << " = "<<ret(c_, i-n) << ", ";
		}
		//Rcout << std::endl;
	}

	X.transpose();

	return ret;
}

NumericMatrix rollAreaSum (Eigen::SparseMatrix<double>, std::vector<double>, int funcID = 1 ,int size = 1e+6 );

//' The numbers start at the first row and end at the last row having a full sized widow
//' 
//' @name rollAreaSum
//' @aliases rollAreaSum,FastWilcoxTest-method
//' @rdname rollAreaSum-methods
//' @docType methods
//' @description calculate a rolling sum of the rows
//' @param X the sparse matrix
//' @param size the size of the rolling window
//' @param location the location for every row in the matrix
//' @param funcID two functions : 1 == sum; 2 == mean
//' @title rolling sum over sparse matrix
//' @export
// [[Rcpp::export]]
NumericMatrix  rollAreaSum (Eigen::SparseMatrix<double> X, std::vector<double> location, 
	int funcID, int size ){

	if ( X.innerSize() != location.size() )
		::Rf_error("Sorry, I need a location for every row of the matrix (m) (%d, %d)", location.size(), X.innerSize() );
	
	//Rcout << "calculating " <<  X.outerSize() << " tests (columns) using "<< CMP.size()<< " resp. " << X.innerSize() << " values" << std::endl;

	std::vector<double> A(X.innerSize());
	std::vector<int> N(X.innerSize());
	int max = 0;
	for ( int i = 0; i < location.size(); i ++){
		for ( int a = 0; a +i < location.size(); a++ ){
			if ( location[a+i] - location[i] >= size){
				N[i] = a;
				max = i;
				//Rcout << "N[" << i<<"] = "<< N[i] <<"; "<< location[a+i] << "-"<< location[i] <<" = "<< (location[a+i] - location[i]) << std::endl;
				break;
			}
		}
	}
	//Rcout << "max=" << max << std::endl;
	NumericMatrix ret( max+1, X.outerSize() );
	std::fill( ret.begin(), ret.end(), 0.0 ) ;

	//Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;
	if ( funcID == 1 ) {
	for (int c_=0; c_ < X.outerSize(); ++c_){
		// copy data from sparse to dense vector
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		//now iterate over all data
		for ( int i=0; i <= max; i++){
			if ( N[i] != 0 ){
				ret(i, c_) = std::accumulate(A.begin()+(i), A.begin()+i+N[i], 0.0);
			}
		}
	}
	}
	else if ( funcID == 2 ) {
		for (int c_=0; c_ < X.outerSize(); ++c_){
			// copy data from sparse to dense vector
			std::fill(A.begin(), A.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
				A[it.row()] =  it.value();
			}
			//now iterate over all data
			for ( int i=0; i <= max; i++){
				if ( N[i] != 0 ){
					ret(i, c_) = std::accumulate(A.begin()+(i), A.begin()+i+N[i], 0.0) / N[i];
				}
			}
		}
	}

	return ret;
}
