// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
#include "FastWilcoxTest_RcppExports.h"


using namespace Rcpp;
#include <R.h>
#include <Rdefines.h>
#include <numeric>

// [[Rcpp::interfaces(r, cpp)]]

//' @name correlationCoefficient
//' @aliases correlationCoefficient,FastWilcoxTest-method
//' @rdname correlationCoefficient-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y
//' @param X one numeric vector
//' @param Y the correlating vector
//' @title Calculate correlation over two double vectors
//' @export
// [[Rcpp::export]]
float correlationCoefficient( std::vector<double> X,  std::vector<double> Y)
{

    double sum_X = 0.0, sum_Y = 0.0, sum_XY = 0.0;
    double squareSum_X = 0.0, squareSum_Y = 0.0;

    if ( X.size() != Y.size() )
    	::Rf_error("Sorry, I need arrays of the same size X and Y(%d, %d)", X.size(), Y.size() );

    int n = X.size();

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X += X[i];

        // sum of elements of array Y.
        sum_Y +=  Y[i];

        // sum of X[i] * Y[i].
        sum_XY += X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X +=  X[i] * X[i];
        squareSum_Y +=  Y[i] * Y[i];
    }

    // use formula for calculating correlation coefficient.
    double corr = (n * sum_XY - sum_X * sum_Y)
                  / sqrt((n * squareSum_X - sum_X * sum_X)
                      * (n * squareSum_Y - sum_Y * sum_Y));

    return corr;
}

// [[Rcpp::interfaces(r, cpp)]]

//' @name CorMatrixIDS
//' @aliases CorMatrixIDS,FastWilcoxTest-method
//' @rdname CorMatrixIDS-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y (slower than apply cor)
//' @param X the sparse matrix
//' @param CMP the vector to correlate every column of the matrix to
//' @param ids the rows of the matrix to correlate to
//' @title Calculate correlation over two double vectors
//' @export
// [[Rcpp::export]]
std::vector<double> CorMatrixIDS (Eigen::MappedSparseMatrix<double> X, std::vector<double> CMP, std::vector<int> ids ) {
	if ( ids.size() != CMP.size() )
		::Rf_error("Sorry, I need arrays of the same size nrow(X[ids,]) and length(CMP)" );

	std::vector<double> A(ids.size());
	std::vector<double> ret(X.cols(), 100);

	ids = FastWilcoxTest::minusOne(ids);
	Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;

	for ( int c_=0; c_ < X.cols(); ++c_ ){
		for ( unsigned int i = 0; i< ids.size(); i++ ) {
			A.at(i) = X.coeff(ids.at(i),c_);
		}
		ret.at(c_) = correlationCoefficient ( CMP, A );
	}
	return ret;
}

//' @name CorMatrix
//' @aliases CorMatrix,FastWilcoxTest-method
//' @rdname CorMatrix-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y
//' approximately 3x faster than an apply using the R cor function on sparse data
//' @param X the sparse matrix
//' @param CMP the vector to correlate every column of the matrix to
//' @title Calculate correlation over two double vectors
//' @export
// [[Rcpp::export]]
std::vector<double>  CorMatrix (Eigen::SparseMatrix<double> X, std::vector<double> CMP){

	X= X.transpose();

	//Rcout << "calculating " <<  X.outerSize() << " tests (columns) using "<< CMP.size()<< " resp. " << X.innerSize() << " values" << std::endl;

	if ( X.innerSize() != CMP.size() )
		::Rf_error("Sorry, I need arrays of the same size ncol(X) and length(CMP) (%d, %d)", X.innerSize(), CMP.size() );

	std::vector<double> A(CMP.size());
	std::vector<double> ret(X.cols());

	//Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;
	for (int c_=0; c_ < X.outerSize(); ++c_){
		// https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		ret.at(c_) = correlationCoefficient ( CMP, A );
	}

	X.transpose();

	return ret;
}


// [[Rcpp::interfaces(r, cpp)]]

//' @name CorNormalMatrix
//' @aliases CorNormalMatrixIDS,FastWilcoxTest-method
//' @rdname CorNormalMatrixIDS-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y 
//' @param X the normal matrix
//' @param CMP the vector to correlate every column of the matrix to
//' @title Calculate correlation over two double vectors
//' @export
// [[Rcpp::export]]
std::vector<double> CorNormalMatrix (NumericMatrix X, std::vector<double> CMP ) {
	if ( X.nrow()  != CMP.size() )
		::Rf_error("Sorry, I need arrays of the same size nrow(X) and length(CMP)(%d, %d)", X.nrow(), CMP.size() );

	std::vector<double> A(X.nrow());
	std::vector<double> ret(X.ncol(), 100);

  	int nr = X.nrow();
	for ( int i =0; i < X.ncol(); i++ ){
		ret.at(i) = correlationCoefficient ( CMP,  std::vector<double>(X.begin() + (nr * i), X.begin() + ((nr * (i+1))) ) );
	}
	return ret;
}


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
