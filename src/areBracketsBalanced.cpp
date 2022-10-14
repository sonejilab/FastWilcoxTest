
#include <Rcpp.h>
//using namespace Rcpp;
#include <iostream>
#include <string>
#include <array>
#include <stdexcept>

//' @name areBracketsBalanced
//' @aliases areBracketsBalanced,FastWilcoxTest-method
//' @rdname areBracketsBalanced-methods
//' @docType methods
//' @description simply check if all brackets are opened and closed in the correct way.
//' A maximum of 500 brackets can be processed by this function.
//' @param xprs the string to check
//' @title check that brackets are used as should
//' @export
//[[Rcpp::export]]
bool areBracketsBalanced( std::string xprs ){

	std::array<int, 3> my_map ={ 0,0,0 };
	std::array<char, 500> brackets;

	int posInB = 0;

	for (int i = 0; i < xprs.length(); i++) {
		if ( posInB ==500 ){
			throw std::out_of_range ("too many brackets -> max 500 !");
		}
		switch ( xprs[i]){
			case '(':
				brackets[posInB++] = xprs[i];
				my_map[0] ++;
				//Rcout << xprs[i] << my_map[0] << std::endl;
				break;
			case ')':
				brackets[posInB++] = xprs[i];
				my_map[0] --;
				//Rcout << xprs[i] << my_map[0] <<std::endl;
				break;
			case '[':
				brackets[posInB++] = xprs[i];
				my_map[1] ++;
				//Rcout << xprs[i] << my_map[1] << std::endl;
				break;
			case ']':
				brackets[posInB++] = xprs[i];
				my_map[1] --;
				//Rcout << xprs[i] << my_map[1] <<std::endl;
				break;
			case '{':
				brackets[posInB++] = xprs[i];
				my_map[2] ++;
				//Rcout << xprs[i] << my_map[2] << std::endl;
				break;
			case '}':
				brackets[posInB++] = xprs[i];
				my_map[2] --;
				//Rcout << xprs[i] << my_map[2] << std::endl;
				break;
		};
	}
	
	for ( int i = 0; i < 3; i++){
		if ( my_map[i] != 0){
			// not all closing are opened or the other way around
			return false;
		}
	}

	for (int i = 0; i < posInB; i++) {
		switch ( brackets[i]){
			case '(':
				if ( brackets[i+1] == ']' || brackets[i+1] == '}'){
					return false;
				}
				break;
			case '[':
				if ( brackets[i+1] == ')' || brackets[i+1] == '}'){
					return false;
				}
				break;
			case '{':
				if ( brackets[i+1] == ')' || brackets[i+1] == ']'){
					return false;
				}
				break;
		};
	}

	return true;

}
