// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
#include <math.h>
using namespace Rcpp;
#include <vector>
#include <stdexcept>
typedef Eigen::MappedSparseMatrix<double> MSpMat;

// [[Rcpp::export]]
double logFC ( std::vector<double> A, std::vector<double> B ) {
	double res = 0.0;
	double Asum = A[0] * 1.0;
	double Bsum = B[0] * 1.0;
	for ( int i=1; i<A.size(); i++ ){
		Asum = log( exp( Asum - A[i]) + 1.0)  + A[i];
	}
	for ( int i=1; i<B.size(); i++ ){
		Bsum = log( exp( Bsum - B[i]) + 1.0 ) + B[i];
	}
	/*Rcout << "Int values A; a size; B; b size:" << Asum <<";"<< A.size()<<";"<< Bsum <<";"<< B.size() << std::endl;*/
	return (Asum - log(A.size()))-(Bsum - log(B.size())) ;
}

std::vector<int> minusOne ( std::vector<int>  X ){
	for ( int i = 0; i < X.size(); i ++) {
		X[i] --;
	}
	return X;
}

// [[Rcpp::export]]
std::vector<double> StatTest (Eigen::MappedSparseMatrix<double> X, std::vector<int> test,
		std::vector<int> backgound, double logFCcut = 1.0, bool display_progress=true ){

    Rcout << "Standard looping over a sparse matrix" << std::endl;

    std::vector<double> logFCpass(X.rows(), 0.0);
    std::vector<double> A(test.size(), 0.0);
    std::vector<double> B(backgound.size(), 0.0 );
    int pass = 0;
    std::vector<int> itA = minusOne( test );
    std::vector<int> itB = minusOne( backgound );

    for ( int c_=0; c_ < X.cols(); c_++ ){
    	for (unsigned int i = 0; i< itA.size(); i++ ) {
    		if ( itA[i] < 0 || itA[i] >= X.rows() ) {
    			 throw std::invalid_argument( "test out of bounds" );
    		}
    		A[i] = X.coeff(itA[i],c_);
    	}
    	for (unsigned int i = 0; i< itB.size(); i++ ) {
    		if (itB[i] < 0 || itB[i] >= X.rows() ) {
    		    			 throw std::invalid_argument( "itB out of bounds" );
    		    		}
    	    B[i] = X.coeff(itB[i],c_);
    	}
    	logFCpass[c_] = logFC( A, B ) > logFCcut;
    	if ( logFCpass[c_] > logFCcut ) {
    		pass++;
    	}
    }

	Rcout << "n return values: " << pass <<std::endl;
	return logFCpass;
}


