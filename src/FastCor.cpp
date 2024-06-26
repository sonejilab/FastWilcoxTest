// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
#include "FastWilcoxTest_RcppExports.h"


using namespace Rcpp;
#include <R.h>
#include <Rdefines.h>
#include <numeric>

#include <progress.hpp>

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
    	::Rf_error("Sorry, I need arrays of the same size X and Y(%lu, %lu)", X.size(), Y.size() );

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
		ret.at(c_) = FastWilcoxTest::correlationCoefficient ( CMP, A );
	}
	return ret;
}

// [[Rcpp::interfaces(r, cpp)]]

//' @name CorMatrixIDS_N
//' @aliases CorMatrixIDS_N,FastWilcoxTest-method
//' @rdname CorMatrixIDS_N-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y (slower than apply cor)
//' @param X the sparse matrix
//' @param CMP the vector to correlate every column of the matrix to
//' @param ids the rows of the matrix to correlate to
//' @title Calculate correlation over two double vectors
//' @returns a matrix of Rho and n (cells where both values where > 0)
//' @export
// [[Rcpp::export]]
NumericMatrix CorMatrixIDS_N (Eigen::MappedSparseMatrix<double> X, std::vector<double> CMP, std::vector<int> ids ) {
	if ( ids.size() != CMP.size() )
		::Rf_error("Sorry, I need arrays of the same size nrow(X[ids,]) and length(CMP)" );

	std::vector<double> A(ids.size());
	NumericMatrix ret(X.cols(), 2);
	
	double zero = 0.0;
	std::fill( ret.begin(), ret.end(), zero ) ;

	ids = FastWilcoxTest::minusOne(ids);
	Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;

	for ( int c_=0; c_ < X.cols(); ++c_ ){
		double n = 0.0;
		for ( unsigned int i = 0; i< ids.size(); i++ ) {
			A.at(i) = X.coeff(ids.at(i),c_);
			n++;
		}
		ret(c_, 0) = FastWilcoxTest::correlationCoefficient ( CMP, A );
		ret(c_, 1) = n;
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
		::Rf_error("Sorry, I need arrays of the same size ncol(X) and length(CMP) (%lu, %lu)", X.innerSize(), CMP.size() );

	std::vector<double> A(CMP.size());
	std::vector<double> ret(X.cols());

	//Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;
	for (int c_=0; c_ < X.outerSize(); ++c_){
		// https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		ret.at(c_) = FastWilcoxTest::correlationCoefficient ( CMP, A );
	}

	X.transpose();

	return ret;
}

//' @name CorMatrix_N
//' @aliases CorMatrix_N,FastWilcoxTest-method
//' @rdname CorMatrix_N-methods
//' @docType methods
//' @description simply calculate the correlation between X and Y
//' approximately 3x faster than an apply using the R cor function on sparse data
//' @param X the sparse matrix
//' @param CMP the vector to correlate every column of the matrix to
//' @title Calculate correlation over two double vectors
//' @returns a matrix of Rho and n (cells where both values where > 0)
//' @export
// [[Rcpp::export]]
NumericMatrix  CorMatrix_N (Eigen::SparseMatrix<double> X, std::vector<double> CMP){

	X= X.transpose();

	//Rcout << "calculating " <<  X.outerSize() << " tests (columns) using "<< CMP.size()<< " resp. " << X.innerSize() << " values" << std::endl;

	if ( X.innerSize() != CMP.size() )
		::Rf_error("Sorry, I need arrays of the same size ncol(X) and length(CMP) (%lu, %lu)", X.innerSize(), CMP.size() );
	std::vector<double> A(CMP.size());
	NumericMatrix ret(X.cols(), 3);
	
	double zero = 0.0;
	std::fill( ret.begin(), ret.end(), zero ) ;

	//Rcout << "calculating " <<  X.cols() << " correlations" << std::endl;
	for (int c_=0; c_ < X.outerSize(); ++c_){
		// https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
		double n = 0.0;
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
			n++;
			
		}

		//Rcout << "n = " <<  n << std::endl;
		if ( n > 2){
			ret(c_, 0) = FastWilcoxTest::correlationCoefficient ( CMP, A );
			ret(c_, 1) = n;
			ret(c_, 2) = ( ret(c_,0) * sqrt( n - 2.0 )) / sqrt( 1.0- ret(c_,0) * ret(c_,0));
		}else {
			ret(c_, 0) = 0.0;
			ret(c_, 1) = n;
			ret(c_, 2) = 0.0;
		}
		
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
		::Rf_error("Sorry, I need arrays of the same size nrow(X) and length(CMP)(%u, %lu)", X.nrow(), CMP.size() );

	std::vector<double> A(X.nrow());
	std::vector<double> ret(X.ncol(), 100);

  	int nr = X.nrow();
	for ( int i =0; i < X.ncol(); i++ ){
		ret.at(i) =  FastWilcoxTest::correlationCoefficient ( CMP,  std::vector<double>(X.begin() + (nr * i), X.begin() + ((nr * (i+1))) ) );
	}
	return ret;
}


//' Calculatethe rolling sum for a max distance from the start
//' The location of each row is taken from the S (Start) vector
//' In R the colnames need to be set to the input colnames whereas the start positions need to be set to chrXY : S[i] - (S[i]+n)
//' @name rollSumStart
//' @aliases rollSumStart,FastWilcoxTest-method
//' @rdname rollSumStartStart-methods
//' @docType methods
//' @description calculate a rolling sum of the rows
//' @param X the sparse matrix
//' @param n the length of the rolling window
//' @param S the start positions of each X row
//' @title rolling sum over sparse matrix
//' @export
// [[Rcpp::export]]
NumericMatrix  rollSumStart (Eigen::SparseMatrix<double> X, double n, std::vector<double> S ){

	Rcout << "calculating log1p(rollSumStart) for " <<  X.outerSize() << " cells " << std::endl;

	Progress p(X.outerSize(), true);

	//X= X.transpose();

	//Rcout << "calculating " <<  X.outerSize() << " tests (columns) using "<< CMP.size()<< " resp. " << X.innerSize() << " values" << std::endl;

	if ( S.back() < n )
		::Rf_error("Sorry, the total dimension of S is smaller than the window dimension) (%f, %f)", S.back(), n );

	//NumericMatrix ret( nrow, ncol );
	std::vector<double> A (X.innerSize());

	double tmp = 0;

	std::vector<int> addN (0);
	std::vector<int> start (0);
	// calculate that only once!
	for ( int i = 0; i<  S.size()-1; i ++ ) {
		for ( int a = i; a < S.size(); a++) {
			// does the distance exceed the max didstance?
			if ( S[a] - S[i] >= n ){
				addN.push_back( a -i );
				start.push_back(i);
//				Rcout << "S[ ->i<- ] "<< S[i] << " is more than n "<< n <<" bp away from S[-> a <-] "<< S[a] << "addN[i] = " << addN[i]<< std::endl;
				i = i + ( a-i ) /2;
				break;
			}
		}
	}
	NumericMatrix ret( addN.size(), X.outerSize() +2 );
	for ( int i = 0; i < addN.size(); i++ ){
		ret(i,0) = start[i]*1.0;
		ret(i,1) = addN[i]*1.0;
	}

	for (int c_=0; c_ < X.outerSize(); ++c_){
		std::fill(A.begin(), A.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		p.increment();
		//Rcout << "summing " <<  X.innerSize() << " values for cell col " << c_ << std::endl;
		//now iterate over all data
		for ( int i=0; i < addN.size(); i++){
			
			//Rcout << "addN = " << addN << std::endl;
			if ( addN[i] < 5 ){
				ret(i, c_+2) = 0.0;
			}
			else {
				ret(i, c_+2) = log1p(std::accumulate(A.begin()+(start[i]), A.begin()+start[i]+addN[i], 0.0));
				// tmp = std::accumulate(A.begin()+(i), A.begin()+i+addN[i], 0.0);
				// if ( tmp > 0 ) {
				// 	ret(c_, i) = tmp / addN[i];
				// }else {
				// 	ret(c_, i) = 0.0;
				// }
			}
			//Rcout << "result = "<<  ret(c_, i) <<" for cell"<< c_ << " and addN = " << addN[i] << " (i="<<i<<" of "<< addN.size() <<")"<< std::endl;
		}
		//break;
		//Rcout << std::endl;
	}

	//X.transpose();
	// colnames(ret) are the colnames of the input matrix and are unknown to c++

	return ret;
}
