#include <Rcpp.h>
#include<iostream>
#include <sstream>


//' return a vector of hex values for a matrix of RGB values (rows)
//' @aliases rgb2hexS,FastWilcoxTest-method
//' @docType methods
//' @name rgb2hexS
//' @param red int value (0..255)
//' @param green int value (0..255)
//' @param blue int value (0..255)
//' @param with_head (bool start with '#' ?)
//' @title hex color RGB values (int)
//' @return string vector of nrow(Mat) length
//' @export
// [[Rcpp::export]]
std::string rgb2hexS(int red, int green, int blue, bool with_head)
{
  std::stringstream ss;
  if (with_head)
    ss << "#";
  ss << std::hex << (red << 16 | green << 8 | blue );
  return ss.str();
}


//' return a vector of hex values for a matrix of RGB values (rows)
//' @aliases rgb2hex,FastWilcoxTest-method
//' @name rgb2hex
//' @docType methods
//' @param Mat the matrix or RGB values (numeric)
//' @title hex color vector 4 RGB matrix
//' @return hex color
//' @export
// [[Rcpp::export]]
Rcpp::StringVector rgb2hex( Rcpp::NumericMatrix Mat ) {
 
   Rcpp::StringVector myvector(Mat.nrow());
   
   for (int i= 0; i< Mat.nrow() ; i++){
      myvector[i] = rgb2hexS( int(Mat(i,0)) , int(Mat(i,1)), int(Mat(i,2)), true);
   }

   return myvector;
}

