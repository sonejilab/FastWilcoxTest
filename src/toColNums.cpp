// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>


using namespace Rcpp;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
std::vector<double> toColNums(Eigen::SparseMatrix<double> data) {
  std::vector<double> tmp( data.nonZeros() );
  int iter=0;
  int id = 0;
  for ( int k=0; k < data.outerSize(); ++k){
	  id ++;
	  for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
	  	  tmp[iter++] = id;
	  }
  }
  return tmp;
}


//' @title ColNotZero returns the amount of not zero values in each column
//' @aliases ColNotZero,FastWilcoxTest-method
//' @rdname ColNotZero
//' @description a c++ implementation of apply(x,2,function(d) {length(which(d!=0))} )
//' @param data a sparse matrix
//' @return a vector with nGene information
//' @export
// [[Rcpp::export]]
std::vector<double> ColNotZero(Eigen::SparseMatrix<double> data) {
  std::vector<double> tmp( data.outerSize() );
  int iter=0;
  int id = 0;
  for ( int k=0; k < data.outerSize(); ++k){
    tmp[iter] = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        tmp[iter] ++;
    }
    iter++;
  }
  return tmp;
}