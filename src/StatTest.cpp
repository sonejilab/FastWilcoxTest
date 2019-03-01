// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
#include <math.h>
#include <stat_rank.h>

using namespace Rcpp;
#include <vector>
#include <stdexcept>
typedef Eigen::MappedSparseMatrix<double> MSpMat;

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#define ABSLOG(x) fabs(log10( (x) ))


// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
double logFC ( std::vector<double> A, std::vector<double> B ) {
	double Asum = A[0] * 1.0;
	double Bsum = B[0] * 1.0;
	for ( unsigned int i=1; i<A.size(); i++ ){
		Asum = log( exp( Asum - A[i]) + 1.0)  + A[i];
	}
	for ( unsigned int i=1; i<B.size(); i++ ){
		Bsum = log( exp( Bsum - B[i]) + 1.0 ) + B[i];
	}
	/*Rcout << "Int values A; a size; B; b size:" << Asum <<";"<< A.size()<<";"<< Bsum <<";"<< B.size() << std::endl;*/
	return (Asum - log(A.size()))-(Bsum - log(B.size())) ;
}

std::vector<int> minusOne ( std::vector<int>  X ){
	for ( unsigned int i = 0; i < X.size(); i ++) {
		X[i] --;
	}
	return X;
}

std::vector<int> plusOne ( std::vector<int>  X ){
	for ( unsigned int i = 0; i < X.size(); i ++) {
		X[i] ++;
	}
	return X;
}

/* direct copy from BioQC/src/wmw_test.c */
/*	greater=0,
	less=1,
	twoSided=2,
	U=3,
	abslog10greater=4,
	log10less=5,
	abslog10twoSided=6,
	Q=7
 */
double wmw_test_stat(double rankSum, int nInds, int nTotal, double tieCoef, int type) {

	double uStat, mu, sigma2, zval, pgt, plt;
	double res;
	int nBg = nTotal-nInds;

	uStat = nInds*nBg+nInds*(nInds+1.0)*0.5-rankSum;

	if(type == 3) {
		res = uStat;
	} else {
		mu = (double)nInds*nBg*0.5; // NOT mu=n1*n2*0.5
		sigma2 = nInds*nBg*(nTotal+1.0)/12.0*tieCoef; //NOT sigma2 = n1*n2*(n+1)/12*tieCoef

		if(type == 0 || type == 4) { /* greater */
			zval = (uStat+0.5-mu)/sqrt(sigma2); // z lower tail
			R::pnorm_both(zval, &pgt, &plt, 0, 0);
			res = type==0 ? pgt : ABSLOG(pgt);
		} else if (type == 1 || type == 5) { /* less */
			zval = (uStat-0.5-mu)/sqrt(sigma2); // z higher tail
			R::pnorm_both(zval, &pgt, &plt, 1, 0);
			res = type==1 ? plt : log10(plt);
		} else if (type == 2 || type == 6 || type == 7) { /* two sided*/
			zval = (uStat-mu-(uStat>mu ? 0.5 : -0.5))/sqrt(sigma2);
			R::pnorm_both(zval, &pgt, &plt, 2, 0);
			res = mu==0.0 ? 1.0 : 2.0*MIN(pgt, plt);
			if(type == 4) {
				res = ABSLOG(res);
			} else if (type == 7) {
				res = pgt<=plt ? ABSLOG(res) : -ABSLOG(res);
			}
		} else {
			/* error("Unrecognized type %d. Should not happen\nPossible only int values  0=greater, 1=less, 2=twoSided, 3=U, 4=abslog10greater, 5=log10less, 6=abslog10twoSided, 7=Q",
            type); */
		}
	}
	return(res);
}

//' @title StatTest runs wilcox test on the columns of the sparse matrix
//' @aliases StatTest,FastWilcoxTest-method
//' @rdname StatTest
//' @description This test implements the Seurat FindMarkers( test.use == "wilcox" ) function
//' in the greatest possible way, but using Rcpp instead of R. So far I could get a ~10x speed improvement.
//' @param X the sparse matrix (tests are applied to columns!)
//' @param interest row IDs for the group of interest
//' @param background row IDS for the background
//' @param logFCcut data is meant to be log() transformed and only columns passing a logFCcut of (default 1) are tested
//' @param minPct only test genes that are detected in a minimum fraction of
//' min.pct cells in either of the two populations. Meant to speed up the function
//' by not testing genes that are very infrequently expressed. Default is 0.1
//' @return a matrix with tested column ids, logFC and p.value
//' @export
// [[Rcpp::export]]
 SEXP StatTest (Eigen::MappedSparseMatrix<double> X, std::vector<int> interest,
		std::vector<int> background, double logFCcut = 1.0, double minPct = 0.1 ){

	//Rcout << "Standard looping over a sparse matrix initializing" << std::endl;
    if ( interest.size() == 0 ){
    	::Rf_error("No values in interest group" );
    }
    if ( background.size() == 0 ){
        ::Rf_error("No values in background group" );
    }
	// internal measurements
	std::vector<double> logFCpass(X.cols(), 0.0);
	std::vector<double> fracInA(X.cols(), 0.0);
	std::vector<double> fracInB(X.cols(), 0.0);
	std::vector<double> indRankSum(X.cols(), 0.0);
    // tmp data storage
	std::vector<double> A(interest.size(), 0.0);
	std::vector<double> B(background.size(), 0.0 );
	double inA = 0;
	double inB = 0;
	double tmp = 0;
    // how many genes pass all filters
	int pass = 0;
	// the corrected (R vs c++) ids
	std::vector<int> itA = minusOne( interest );
	std::vector<int> itB = minusOne( background );

	Rcout << "calculating filters logFC and minPct" << std::endl;
	for ( int c_=0; c_ < X.cols(); ++c_ ){
		inA = 0;
		inB = 0;
		for ( unsigned int i = 0; i< itA.size(); i++ ) {
			if ( itA.at(i) < 0 || itA.at(i) >= X.rows() ) {
				::Rf_error( "test out of bounds" );
			}
			tmp = X.coeff(itA.at(i),c_);
			if ( tmp > 0 ){
				inA ++;
			}
			A[i] = tmp;
		}
		for ( unsigned int i = 0; i< itB.size(); i++ ) {
			if (itB.at(i) < 0 || itB.at(i) >= X.rows() ) {
				::Rf_error("itB out of bounds" );
			}
			tmp = X.coeff(itB.at(i),c_);
			if ( tmp > 0 ){
				inB ++;
			}
			B.at(i) = tmp;
		}
		fracInA[c_] = inA / interest.size();
		fracInB[c_] = inB / background.size();
		logFCpass[c_] = logFC( A, B );
		if ( logFCpass.at(c_) > logFCcut && ( (fracInA[c_] > minPct) +  (fracInB[c_] > minPct) ) > 0 ) {
			pass++;
		}
	}
	if ( pass == 0 ){
			::Rf_error("No gene passed the logFC + min expressed filter - try changing the minPct and logFCcut variables" );
	}
	else {
		Rcout << "calculating wilcox test(s) for " << pass << " genes" << std::endl;
	}

	/* allocate a result 'matrix' */
	NumericMatrix res(pass, 6);
	int n = X.rows();
	std::vector<double> total( itA.size() + itB.size() , -1.0 );

	int id = 0;
	int j;
	int nInd;
	double tie;
	//double indRankSum;

	DRankList list;

	for ( int c_=0; c_ < X.cols(); c_++ ){
		if ( logFCpass[c_] > logFCcut ) {

			/*Test stats copied from the BioOC package */
			j = 0;
			for (unsigned int i = 0; i< itA.size(); i++ ) {
				total.at(j++) = X.coeff(itA.at(i) ,c_);
			}
			for (unsigned int i = 0; i< itB.size(); i++ ) {
				total.at(j++) = X.coeff(itB.at(i) ,c_);
			}
			n = j;

			// populate the DRankList object

			list.refill(total, n);
			//Rcout << "isRanked " << list.isRanked() <<std::endl;
			//list.sortRankDRankList();
			//Rcout << "isRanked " << list.isRanked() <<std::endl;
			// test DRankLint internals
			//continue;

			list.prepareDRankList();

			//Rcout << "finished with the prepare:" << std::endl;
			//list.print();

			tie = list.tieCoef;

			nInd=itA.size();

			indRankSum.at(c_) = 0.0;
			//Rcout << "calculating for gene id " << c_ << " indRankSum starting at " << indRankSum.at(c_) << std::endl;
			for(j=0; j<nInd; ++j) {
				if(!(itA.at(j)>=0 && itA.at(j)<=n-1))
					::Rf_error("Index out of range: gene set %d, gene %d\n", c_+1, j+1);
				if ( list.list.size() <= itA.at(j) )
					::Rf_error("Not enough values in the ranked list; list size %d <= data size", list.list.size(), itA.at(j));
				//Rcout << "adding to RankSum for gene id " << c_ << " and index == "
				//		<<  itA.at(j) << " and value " <<  list.list.at(itA.at(j)).rank << std::endl;
				indRankSum.at(c_) += list.list.at(itA.at(j)).rank;
			}
			// never destroy the  DRankList - that kills the R gc() functionality!!
			//Rcout << "got a result for id " << c_ << " indRankSum == " << indRankSum.at(c_) << std::endl;
			/* store the results */
			res(id,0) = c_ + 1;
			res(id,1) = logFCpass.at(c_);
			res(id,2) = fracInA.at(c_);
			res(id,3) = fracInB.at(c_);
			res(id,4) = indRankSum.at(c_);
			/* store the higher p value as we do drop all lower anyhow. */
			res(id,5) = wmw_test_stat(indRankSum.at(c_), nInd, n, tie, 0);
			//Rcout << "got a result for id " << c_ << "p.value == " << res(id,5) << std::endl;
			id ++;
		}
	}
	colnames(res) = CharacterVector::create("colID", "logFC", "fracExprIN", "fracExprOUT", "rank.sum", "p.value");
	Rcout << "n return values: " << pass <<std::endl;
	return res;
}


