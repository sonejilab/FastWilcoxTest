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
#include <progress.hpp>
#include "FastWilcoxTest_RcppExports.h"
#include <FastCor.h>

// [[Rcpp::interfaces(r, cpp)]]

//' @name LinLang
//' @aliases LinLang,FastWilcoxTest-method
//' @rdname LinLang-methods
//' @docType methods
//' @description Identify genes slowly rising above the detection limit
//' @param X the sparse Matrix (row = genes, col = cells)
//' @param Grouping a numeric vector of group IDs
//' @param nGroup the number of groups
//' @param display_progress show a progress bar (TRUE)
//' @title LinLang test for rise above the detection limit
//' @export
//[[Rcpp::export]]
std::vector<double> LinLang(Eigen::SparseMatrix<double> X, std::vector<int> Grouping, int nGroup, bool display_progress=true ){

	Grouping = FastWilcoxTest::minusOne(Grouping);
	std::vector<double> A(X.innerSize());
	std::vector<double> ret(X.outerSize());

	std::vector<double> means(nGroup);
	std::vector<double> amount(nGroup);
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
			if ( A[i_] != 0 ) {
				amount[Grouping[i_]] ++;
				means[Grouping[i_]] += A[i_];
			}
		}
		for ( int i_ = 0; i_ < nGroup; i_ ++){
			if ( amount[i_] > 0 ){
				means[i_] = means[i_] / amount[i_];
				amount[i_] = amount[i_] / total[i_];
			}

		}
		// now calculate the correlation value and be done with that!
		ret[c_] = FastWilcoxTest::correlationCoefficient( amount, means );
	}
	return(ret);
}