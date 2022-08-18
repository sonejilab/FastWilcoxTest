// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <Rdefines.h>

using namespace Rcpp;

typedef Eigen::Triplet<double> T;

//' @title reshuffle data based on a sparse matrix assuming max double the amount of entries not being zero
//' @aliases ShuffleMatrix,FastWilcoxTest-method
//' @rdname ShuffleMatrix
//' @description replacing the synthetic1 function of RFclust.SGE package 
//' @param X the sparse matrix (tests are applied to columns!)
//' @param maxCols the amount of random columns to send back (default 50)
//' @return a matrix with x, j and i avalues to be put into a new sparse matrix
//' @export
// [[Rcpp::export]]
 Eigen::SparseMatrix<double> ShuffleMatrix (Eigen::SparseMatrix<double> X, int maxCols = 50) {

 	std::vector<T> trp;
 	int max = X.nonZeros();
 	int maxNewGene = 0;
 	trp.reserve( X.nonZeros() );
 	if ( X.outerSize() < maxCols){
 		maxCols = X.outerSize();
 	}

	//std::vector<double> A(X.innerSize());
	NumericVector A(X.innerSize());
	std::vector<bool> geneTaken(X.innerSize());

	//NumericVector Resampled(1);
	int randGene = 0;
	int prob =1;
	int i;
	int trpID = 0;
	int tmp = 0;
	int goodA = 0;
	//Rcout << "create random triplets n=" << max << std::endl;

	for (int c_=0; c_ < maxCols; ++c_){
		std::fill(A.begin(), A.end(), 0.0);
		std::fill(geneTaken.begin(), geneTaken.end(), false);

		i = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[i++] =  it.value();
			//Rcout << "old value " <<  A[i-1] << std::endl;
		}
		goodA= i;
		//Rcout << "calculate randoms " <<  A.length() << std::endl;
		//NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
		// and now we need to fill that into the sparse matrix..
		
		for ( int i=0; i < goodA; i++){
			if ( A[i] == 0.0 || A[i] == 0){
				//Rcout << "col broken at A["<<i<<"]="<<A[i]<< " and running id= " <<  trpID << std::endl;
				break;
			}

			randGene = std::floor(R::runif(0,1) * A.length()) ;

			while ( geneTaken[randGene] ){
				randGene = std::floor(R::runif(0,1) * A.length()) ;
			}
			geneTaken[randGene] = true;

			if ( i == 0 ){
				//Rcout << "inst " <<  randGene << ","<< c_ <<" value " << A[i] << std::endl;
				tmp = trpID;
			}
			
			if ( maxNewGene < randGene){
				maxNewGene = randGene ;
			}
			if ( trpID > max ){
				//Rcout << "triplets exceed max counts - stop" << std::endl;
			}else {
				trp.push_back(T(randGene, c_, A[i]) ) ;
				trpID++;
			}
		}
		//Rcout << "max geneID = "<< maxNewGene << std::endl;
		//Rcout << "do we have access to the triplets? t[" << tmp << "] = " << trp[tmp].row() <<
		//		", " <<  trp[tmp].col() << ", "<< trp[tmp].value() << std::endl;
	}

	//Rcout << "inserted a total of " << trpID << " values into the triplet list"<< std::endl;
	Eigen::SparseMatrix<double> mat2(X.innerSize(), maxCols );

	mat2.setFromTriplets(trp.begin(), trp.end());
	//Rcout << "inserted a total of " << inserted << " values into the dgCMatrix"<< std::endl;
	return ( mat2 );
}