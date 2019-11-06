#include <R.h>
//#include <Rmath.h>
//#include <Rdefines.h>
#include <vector>
#include <stat_rank.h>
#include <cmath>


// [[Rcpp::interfaces(r, cpp)]]
//' Calculate the euclidian distance between consecutive points
//' Can also produce the integral values of the distance.
//' @name euclidian_distances
//' @aliases euclidian_distances,FastWilcoxTest-method
//' @rdname euclidian_distances-methods
//' @docType methods
//' @description calculates the (2D) euclidian distance for a set of x and y values
//' @param X one ORDERED numeric vector
//' @param Y the other vector
//' @param sum create a total sum of these values (integral) default = FALSE
//' @title Calculate  over two double vectors
//' @export
// [[Rcpp::export]]
std::vector<double> euclidian_distances( std::vector<double> X,  std::vector<double> Y, bool sum = false)
{
	std::vector<double> distance ( X.size() );
	distance[0] = 0;
	for ( int i = 1; i< X.size(); i++ ) {
		distance[i] = sqrt(	pow( (X[i-1] - X[i]) ,2 ) +	pow( (Y[i-1] - Y[i]) ,2 )  );
	}
	if ( sum ) {
		for ( int i = 1; i< X.size(); i++ ) {
			distance[i] += distance[i-1];
		}
	}
	return distance;
}

// [[Rcpp::interfaces(r, cpp)]]
//' Calculate the euclidian distance between consecutive points
//' Can also produce the integral values of the distance.
//' @name euclidian_distances3d
//' @aliases euclidian_distances3d,FastWilcoxTest-method
//' @rdname euclidian_distances3d-methods
//' @docType methods
//' @description calculates the (3D) euclidian distance for a set of x and y values
//' @param X one ORDERED numeric vector
//' @param Y the other vector
//' @param Z the thrid dimension
//' @param sum create a total sum of these values (integral) default = FALSE
//' @title Calculate  over two double vectors
//' @export
// [[Rcpp::export]]
std::vector<double> euclidian_distances3d( std::vector<double> X,  std::vector<double> Y, std::vector<double> Z, bool sum = false)
{
	std::vector<double> distance ( X.size() );
	distance[0] = 0;
	for ( int i = 1; i< X.size(); i++ ) {
		distance[i] = sqrt(	pow( (X[i-1] - X[i]) ,2 ) +	pow( (Y[i-1] - Y[i]) ,2 ) + pow( (Z[i-1] - Z[i]) ,2 )  );
	}
	if ( sum ) {
		for ( int i = 1; i< X.size(); i++ ) {
			distance[i] += distance[i-1];
		}
	}
	return distance;
}


// [[Rcpp::interfaces(r, cpp)]]
//' use the eucledian distance between one cell and all cells to find the order in the data
//' @name eDist3d
//' @aliases eDist3d,FastWilcoxTest-method
//' @rdname eDist3d-methods
//' @docType methods
//' @description calculates the (3D) euclidian distance for a set of x and y values
//' @param X one numeric vector
//' @param Y the other vector
//' @param Z the thrid dimension
//' @param id which id to start from
//' @title find the euclidian order in a 3D vector
//' @export
// [[Rcpp::export]]
std::vector<double> eDist3d( std::vector<double> X,  std::vector<double> Y, std::vector<double> Z, int id) {
	std::vector<double> distance ( X.size() );
	std::fill(distance.begin(), distance.end(), 100000000.0);
	for ( int i = 0; i< X.size(); i++ ) {
		//if ( ! X[i] == 0.0 & Y[i] == 0.0 ) {
			distance[i] = sqrt(	pow( (X[i] - X[id]) ,2 ) +	pow( (Y[i] - Y[id]) ,2 ) + pow( (Z[i] - Z[id]) ,2 )  );
		//}
	}
	return(distance);
}
