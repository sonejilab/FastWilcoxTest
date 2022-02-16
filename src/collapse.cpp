
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppEigen.h>
using namespace Rcpp;
#include "FastWilcoxTest_RcppExports.h"

//' @name collapse
//' @aliases collapse,FastWilcoxTest-method
//' @rdname collapse-methods
//' @docType methods
//' @description sums up the values for each ids type
//' @param X the sparse matrix
//' @param ids group ids (int vector from 1 10 maxgroup for each column)
//' @param type ( 0: logAdd (defunct); 1 : simple addition; 2: mean )
//' @title collapse the data collumns based on the ids info
//' @export
// [[Rcpp::export]]
NumericMatrix collapse (Eigen::SparseMatrix<double> X, std::vector<int> ids, int type ) {
	int total = 0;
	std::vector<int> USE( X.cols() );

	if ( ids.size() != X.cols() )
		::Rf_error("the ids vector needs to have the same length as the matrix columns ", type);
	for ( int i = 0; i < ids.size(); i++ ){
		if ( total < ids[i])
			total = ids[i];
	}
	ids = FastWilcoxTest::minusOne(ids);

	Rcout << "merge a matrix of " << X.rows() << " rows and "<<  X.cols() <<" into " << total << " return columns";
	switch(type) {
		case 0 :
			Rcout << " using logAdd"  << std::endl;
			::Rf_error("Sorry that is not implemented correctly" );
			break;
		case 1 :
			Rcout << " using simple add"  << std::endl;
			break;
		case 2 :
			Rcout << " mean"  << std::endl;
			std::fill( USE.begin(), USE.end(), 0 ) ;
			for ( int i = 0; i < ids.size(); i++ ){
				USE.at(ids[i]) ++;
			}
			break;
		default: ::Rf_error("Sorry the merge function type %d is not defined ", type);
	}
	X= X.transpose();

	NumericMatrix res(X.cols(), total );
	Rcout << " I create a results table with "  <<  X.cols() << "rows and " << total << " columns" << std::endl;
	std::fill( res.begin(), res.end(), 0.0 ) ;
	int tmp;

	for (int c_=0; c_ < X.outerSize(); ++c_){
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			switch(type) {
			case 0 :
				res( c_, ids[it.row()]) += exp(it.value());
				break;
			case 1 :
				res( c_, ids[it.row()]) += it.value();
				break;
			case 2 :
				res( c_, ids[it.row()]) += it.value();
				break;
			}
		}
	}
	int nrows = res.nrow();
	int ncolumns = res.ncol();
	switch(type) {
		case 0 :
			for (int c_=0; c_ < ncolumns; c_++){
				for (int i = 0; i < nrows; i++){
					res( i, c_) = log( res( i, c_ ) +1 );
				}
			} 
			break;
		case 2 :
			for (int c_=0; c_ < ncolumns; c_++){
				for (int i = 0; i < nrows; i++){
					res( i, c_) = res( i, c_ ) / USE[c_];
				}
			}
			break;
	}

	

	return res;
}
