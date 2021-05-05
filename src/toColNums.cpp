// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>
#include <iostream>
#include <fstream> 

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

//' @title sparse2SQLite_text_file creates a simple text file from the matrix contents
//' @aliases sparse2SQLite_text_file,FastWilcoxTest-method
//' @rdname sparse2SQLite_text_file
//' @description circumvent the memory expenses during RSQlite database creation (melting the whole matrix)
//' @param data a sparse matrix
//' @return a vector with nGene information
//' @export
// [[Rcpp::export]]
void sparse2SQLite_text_file( Eigen::SparseMatrix<double> data, String file, char sep=' '){
  std::ofstream o(file);

  for ( int k=0; k < data.outerSize(); ++k){

    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      /*
      gene_id= x@i +1, 
      cell_id = toColNums( x ), 
      value= x@x
      */
      o << it.row() + 1 << sep << it.col() + 1 << sep << it.value() << std::endl;
    }
  }
  o.close();
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