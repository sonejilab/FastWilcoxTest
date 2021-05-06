// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <fstream> 

using namespace Rcpp;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]

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

  //char buf[512000];
  //o.rdbuf()->pubsetbuf(buf,512000);


  for ( int k=0; k < data.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      o << it.row() + 1 << sep << it.col() + 1 << sep << it.value() << "\n";
    }
  }
  o.close();
}
