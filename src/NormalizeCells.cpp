/* Replicate the 'Seurat' 'RunUMISamplingPerCell' function using the sparse nature of the matrix.
 * This function is used in the BioData class (https://github.com/stela2502/BioData)
 *
 */

#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>
using namespace Rcpp;
#include <progress.hpp>
//#include <Rcpp::sugar::sample.h>
#include <stat_rank.h>



//' @name NormalizeCells
//' @aliases NormalizeCells,FastWilcoxTest-method
//' @rdname NormalizeCells-methods
//' @docType methods
//' @description Normalize the single cell expression values to a total of nUMI reads.
//' @param data the sparse Matrix (row = genes, col = cells)
//' @param nUMI aim to normalize each cell to (cells expressing less are set to 0)
//' @param display_progress show a progress bar (TRUE)
//' @title UMI normalize a single cell expression matrix
//' @export
//[[Rcpp::export]]
Eigen::SparseMatrix<double>  NormalizeCells (Eigen::SparseMatrix<double> X, int nUMI, bool display_progress=true ){



	Progress p(X.outerSize(), display_progress);


	double sum;
	double fact; // the scaling factor
	const double zero=0.0; // to set lost cells to 0
	const double l = -1.0; // to set lost genes to -1
	double nK, rn; // the below 0 fraction of the normalized value
	int less_max, more_max;
	DRankList lostAndGained( X.rows() ); // far too big list

	//Rcout << "inizialized  DRankList with " << X.rows() << " possible values" << std::endl;
	nK = 0;
	for (int k=0; k < X.outerSize(); ++k){
		sum = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it){
			sum += it.value();
		}
		if ( nK < sum )
			nK = sum;
	}
	if ( nUMI < (nK / 10) ){
		p.cleanup();
		::Rf_error( "Please adjust the nUMI value - you would unnecessarily lose about 90 percent of the available information" );
	}

	for (int k=0; k < X.outerSize(); ++k){

		p.increment();
		sum = 0;
		less_max = 0; more_max = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it){
			sum += it.value();
		}
		if ( sum == nUMI) {
			continue;
		}
		// Rcout << "sum " << sum;
		if (sum < nUMI ) { // drop this cell (set all to 0)
			// Rcout << " is below " <<  nUMI << " set to 0" << std::endl;
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it){
				it.valueRef() = zero;
			}

		}else {
			fact = double(sum / nUMI);
			sum = 0.0;

			lostAndGained.len=X.rows();
			lostAndGained.purge();
			lostAndGained.ranked = true; //I do not need the ranking function here!

			 //Rcout <<"sum " << sum<< " is OK fact set to " <<  fact << std::endl;
			/* new logics - store the random difference between nK and rn
			* and later on select those with the highest negative/positive difference.
			* To make that work store 2-diff respectively 1-diff for the <0.2 nK resp <0.5 nK that still passed and
			* and -2+diff resp -1+diff for the >8.8 nK respectively >= 0.5 that failed the test.
			*/
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it){
				double neu;
				nK = modf((it.value() /fact), &neu);
				rn = R::runif(0,1);
				//Rcout << "makes sense?: value " <<it.value() << "/"<< fact << " (fact) = " << neu << "." << nK << std::endl;

				if (  rn  <= nK )  { //the test passed
					neu ++; // new value +1
					if ( nK < 0.2 ) {
						//low, but passed - could be dropped if we get too many
						lostAndGained.add( it.row(), ( 3.0 - ( rn - nK ) ) );
						more_max++;
					}else if ( nK <= 0.5 ) {
						lostAndGained.add( it.row(), ( 2.0 - ( rn - nK ) ) );
						more_max++;
					}else {
						lostAndGained.add( it.row(), (  rn - nK  ) ); //positive
						more_max++;
					}
				}else{ //test failed
					if ( nK > 0.8 ) {
						lostAndGained.add( it.row(), ( -3.0 + rn - nK) );
						less_max ++;
					}else if ( nK >= 0.5 ) {
						lostAndGained.add( it.row(), ( -2.0 + rn - nK)  );
						less_max ++;
					}else {
						lostAndGained.add( it.row(), (  nK - rn )  ); //negative
						less_max ++;
					}
				}
				if ( neu <= 0.0 ) {
					it.valueRef() = l; // lost expression (-1.0)
				}else {
					it.valueRef() = neu;
					sum += neu;
				}
			}
			if ( sum < nUMI ) {
				// reduce some values
				int diff = nUMI - sum +1;
				if ( more_max < diff )
					::Rf_error( "Not enough close values to fix the  sum (%d) < nUMI problem (%d) (max $d)" , sum, nUMI, more_max);

				//Rcout << k << " sum is smaller than nUMI - need to add +1 x " << diff << " values (sum=" <<  sum<<")" << std::endl;
				lostAndGained.sortDRankList(); // sort by it value
				//Rcout << "list sorted " <<std::endl;
				for( int i = 0; i < diff; i++){
					//Rcout << "   value before (+) " << X.coeff( lostAndGained.list.at( i ).index ,k ) << " ("<<i<<") ";
					if ( lostAndGained.list.at( i ).index == -1 ){
						diff ++;
						continue;
					}
					if ( X.coeff( lostAndGained.list.at( i ).index ,k ) == -1.0 )
						X.coeffRef( lostAndGained.list.at( i ).index ,k ) = 1.0;
					else
						X.coeffRef( lostAndGained.list.at( i ).index ,k ) = X.coeff( lostAndGained.list.at( i ).index ,k ) +1.0;
					//Rcout << "after " << X.coeff( lostAndGained.list.at( i ).index ,k ) << std::endl;
				}
			}else if ( sum > nUMI ) {
				int diff =  sum -nUMI;
				if ( less_max < diff )
					::Rf_error( "Not enough close values to fix the  sum (%d) > nUMI (%d) problem (max %d)" , sum, nUMI, less_max);

				//Rcout << k <<" sum is bigger than nUMI - need to remove -1 x " << diff << " values (sum=" <<  sum<<")" << std::endl;

				lostAndGained.sortDRankList(); // sort by it value
				for( int i = 0; i < diff; i++){
					//Rcout << "   value before (-) " << X.coeff( lostAndGained.list.at( lostAndGained.len - i ).index ,k ) << " ("<< ( lostAndGained.len - i )<< ") ";
					if ( lostAndGained.list.at( lostAndGained.len - i ).index == -1 ){
						diff ++;
						continue;
					}
					if ( X.coeff( lostAndGained.list.at( lostAndGained.len - i ).index ,k ) == -1.0 ){
						diff ++;
					}else {
						nK = X.coeffRef( lostAndGained.list.at( lostAndGained.len - i ).index ,k );
						nK -= 1.0;
						if ( nK  == 0.0 ){
							X.coeffRef( lostAndGained.list.at( lostAndGained.len - i ).index ,k ) = l; //lost
						}else {
							X.coeffRef( lostAndGained.list.at( lostAndGained.len - i ).index ,k ) = nK;
						}
					}
					//Rcout << "after " << X.coeff( lostAndGained.list.at( lostAndGained.len -i ).index ,k ) << std::endl;
				}
			}
		}
	}
	//p.cleanup();
	return (X);
}
