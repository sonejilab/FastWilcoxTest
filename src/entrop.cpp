
#include <progress.hpp>
#include <string.h>
#include <FastWilcoxTest.h>
#include <math.h> 


// [[Rcpp::interfaces(r, cpp)]]
//' @name entropy
//' @aliases entropy,cellexalvrR-method
//' @rdname entropy-methods
//' @docType methods
//' @description calculate the entropy of one double vector. Each number is a new group.
//' @param X a double vector of group ids
//' @title description of function entrop
//' @export 
// [[Rcpp::export]]
double entropy( std::vector<double> X )
{
	Table data;
	return ( data.Entropy( X ) );
}


/*
re-implement this functionallity in c++!

for( res in  lapply( 1:length(gvect), function( id ) {
                d = FastWilcoxTest::eDist3d( x[,1], x[,2],x[,3], id -1 )
#               order(d)[1:n]
                pb$tick()
                m = matrix(
                   unlist( lapply( n , function( N ) {
                        OK = which(d < N)
                        if ( length(OK) == 0 ){q
                                return(0)
                        }
                        c( FastWilcoxTest::entropy(  gvect[ OK ]), length(OK) )
                   } ))
                , nrow=2)
                rownames(m) = paste( 'cell', id, c('entropy', 'selectedCells'))
                m

        })) {
              closest = rbind( closest, res )
        }
/**/

// [[Rcpp::interfaces(r, cpp)]]
//' @name SphericEntropy
//' @aliases SphericEntropy,cellexalvrR-method
//' @rdname SphericEntropy-methods
//' @docType methods
//' @description Calculate a per cell entropy for a vector of radius cut-offs.
//' @param X1 a double vector of dim1
//' @param X2 a double vector of dim2
//' @param X3 a double vector of dim3
//' @param gvect a double vector of the output group
//' @param radii a double vector of the euclidian max distances (radius) to be checked
//' @title description of function entrop
//' @export 
// [[Rcpp::export]]
NumericMatrix SphericEntropy( std::vector<double> X1, std::vector<double> X2, std::vector<double> X3, std::vector<double> gvect, std::vector<double> radii ) {
	NumericMatrix res (  X1.size() * 2,  radii.size() );
	std::vector<double> d;
	std::vector<Table> Tables ( radii.size() );

	Progress p(X1.size(), true);

	CharacterVector rh( X1.size()*2 );

	for ( unsigned i = 0; i < X1.size(); i++ ){
		d = FastWilcoxTest::eDist3d( X1, X2, X3, i );
		for ( int a = 0; a < radii.size(); a++){
			Tables.at(a).reset();
		}
		for ( int cellID = 0; cellID < d.size(); cellID ++ ){
			for ( int a = radii.size() -1; a >=0; a --) {
				//Rcout << "test on line" << i << "/"<< cellID<< " : " << d.at(i) << " < " <<radii.at(a) ;
				if ( d.at(cellID) < radii.at(a)){
					//Rcout << "adding " << gvect.at(cellID)<< " at level "<< a << std::endl;
					Tables.at(a).add( gvect.at(cellID) );
				}else {
					//Rcout << " - NO -> break!" << std::endl;
					break;
				}
			}
		}
		for ( int a = 0; a < radii.size(); a++){
			//Tables.at(1).print();
			res(i*2,a) = (double)Tables.at(a)._entropy();
			res(i*2+1,a) = (double)Tables.at(a).sum();
			//Rcout << "cell " << i << " entropy " << res(i*2,a) << " n cell "<< res(i*2+1,a) <<std::endl;
			rh.at(i*2) = "cell " + std::to_string(i+1) + " entropy";
			rh.at(i*2+1) = "cell " + std::to_string(i+1) + " selectedCells";
		}
		p.increment();
	}
	//Rcout << "working fine!" << std::endl;
	CharacterVector ch(radii.size());
	for( int i = 0; i < radii.size(); i++){
		ch(i) = std::to_string(radii.at(i));
	}
	colnames(res) = ch;
	rownames(res) = rh;
	return ( res );
}
