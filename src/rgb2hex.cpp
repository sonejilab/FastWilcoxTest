#include <Rcpp.h>
#include<iostream>
#include <sstream>


// https://www.geeksforgeeks.org/convert-the-given-rgb-color-code-to-hex-color-code/

std::string decToHexa(int n)
{
    // char array to store hexadecimal number
    char hexaDeciNum[2];
 
    // counter for hexadecimal number array
    int i = 0;
    while (n != 0) {
 
        // temporary variable to store remainder
        int temp = 0;
 
        // storing remainder in temp variable.
        temp = n % 16;
 
        // check if temp < 10
        if (temp < 10) {
            hexaDeciNum[i] = temp + 48;
            i++;
        }
        else {
            hexaDeciNum[i] = temp + 55;
            i++;
        }
 
        n = n / 16;
    }
 
    std::string hexCode = "";

    if (i == 2) {
        hexCode.push_back(hexaDeciNum[1]);
        hexCode.push_back(hexaDeciNum[0]);
    }
    else if (i == 1) {
        hexCode = "0";
        hexCode.push_back(hexaDeciNum[1]);
    }
    else if (i == 0)
        hexCode = "00";
 
    // Return the equivalent
    // hexadecimal color code
    return hexCode;
}
 


//' return a vector of hex values for a matrix of RGB values (rows)
//' @aliases rgb2hexS,FastWilcoxTest-method
//' @docType methods
//' @name rgb2hexS
//' @param red int value (0..255)
//' @param green int value (0..255)
//' @param blue int value (0..255)
//' @title hex color RGB values (int)
//' @return string vector of nrow(Mat) length
//' @export
// [[Rcpp::export]]
std::string rgb2hexS(int red, int green, int blue)
{

    if ((red >= 0 && red <= 255)
        && (green >= 0 && green <= 255)
        && (blue >= 0 && blue <= 255)) {
 
        std::string hexCode = "#";
        hexCode += decToHexa(red);
        hexCode += decToHexa(green);
        hexCode += decToHexa(blue);
 
        return hexCode;
    }
 
    // The hex color code doesn't exist
    else
        return "-1";
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
      myvector[i] = rgb2hexS( int(Mat(i,0)) , int(Mat(i,1)), int(Mat(i,2)));
   }

   return myvector;
}

