/* This function will check in a bootstrap approach whether the mean expression
 * level of a gene is linearly correlated to the fraction of cell showing expression of the gene.
 *
 * The rational for the test is that if a clustering is good, any given gene, that is expressed differentially
 * in this grouping should show a (steady) rise above the detection limit.
 * Below the detection limit the gene would not be detected, but the probability of detection would rise.
 * The linear correlation of both the fraction of positive cells and the mean expression level should
 * reveal this rise above the detection level.
 */

#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>
using namespace Rcpp;
#include "FastWilcoxTest_RcppExports.h"
#include <FastCor.h>
#include <Rcpp/sugar/functions/functions.h>


double var ( std::vector<double> v ) {
	double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m =  sum / v.size();

	double accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
    	accum += (d - m) * (d - m);
	});
	double var = accum / (v.size()-1);
	return var;
}

// [[Rcpp::interfaces(r, cpp)]]

//' @name LinLang
//' @aliases LinLang,FastWilcoxTest-method
//' @rdname LinLang-methods
//' @docType methods
//' @description Identify genes slowly rising above the detection limit
//' @param X the sparse Matrix (row = genes, col = cells)
//' @param Grouping a numeric vector of group IDs
//' @param nGroup the number of groups
//' @param minPct ignore genes with less than a fraction of 0.1 (default) of the cells expressing them
//' @param display_progress show a progress bar (TRUE)
//' @title LinLang test for rise above the detection limit
//' @export
// [[Rcpp::export]]
NumericMatrix LinLang (Eigen::SparseMatrix<double> X, std::vector<int> Grouping, int nGroup, double minPct = 0.1, bool display_progress=true ){

	Grouping = FastWilcoxTest::minusOne(Grouping);
	std::vector<double> A(X.innerSize());
	//std::vector<double> ret(X.outerSize());
	NumericMatrix ret( X.outerSize(), 5 );
	std::fill(ret.begin(), ret.end(), 0.0);
	
	
    double sum_of_elems;
    
	std::vector<double> means(nGroup);
	std::vector<double> amount(nGroup);
	std::vector<double> cmp(nGroup);
	// fill cmp with 1:nGroup
	std::iota(cmp.begin(), cmp.end(), 1);
	std::vector<double> total(nGroup);

	std::fill(total.begin(), total.end(), 0.0);
	for ( int i_ = 0; i_ < Grouping.size(); i_ ++){
		total[Grouping[i_]] ++;
	}

    Progress p(X.outerSize(), display_progress);

	for (int c_=0; c_ < X.outerSize(); ++c_){
		p.increment();
		// https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
		std::fill(A.begin(), A.end(), 0.0);
		std::fill(means.begin(), means.end(), 0.0);
		std::fill(amount.begin(), amount.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
			A[it.row()] =  it.value();
		}
		// collect variables
		for ( int i_ = 0; i_ < A.size(); i_ ++ ) {
			if ( A[i_] > 0 ) {
				amount[Grouping[i_]] ++;
				means[Grouping[i_]] += A[i_];
			}
		}
		sum_of_elems = std::accumulate(amount.begin(), amount.end(), 0.0);
		//Rcout << "sum_of_elems " << sum_of_elems << std::endl;
		if ( sum_of_elems / X.innerSize()  > minPct ) {
			//Rcout << "Processed!" << std::endl;
			for ( int i_ = 0; i_ < nGroup; i_ ++){
				if ( amount[i_] > 0 ){
					means[i_] = means[i_] / amount[i_];
					amount[i_] = amount[i_] / total[i_];
				}	
			}
			// now calculate the correlation value and be done with that!
			ret(c_, 3) = var(amount);
			ret(c_, 4) = var(means);
			if ( ret(c_, 3) != 0 ){
				ret(c_, 0) = FastWilcoxTest::correlationCoefficient( cmp, amount);
			}
			if ( ret(c_, 4) != 0 ){
				ret(c_, 1) = FastWilcoxTest::correlationCoefficient( cmp,means );
			}
			if ( ret(c_, 3) != 0 & ret(c_, 4) != 0 ) {
				ret(c_, 2) = FastWilcoxTest::correlationCoefficient( amount, means );
			}
		}
	}
	colnames(ret) = CharacterVector::create("OnOff vs order", "Level vs order", "OnOf vs Level", "OnOf var", "Level var" );
	return ret;
}