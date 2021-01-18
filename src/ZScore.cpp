// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>
using namespace Rcpp;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;
using  Eigen::Map;
using  Eigen::VectorXd;
typedef  Map<VectorXd>  MapVecd;

//' @name ZScore
//' @aliases ZScore,FastWilcoxTest-method
//' @rdname ZScore-methods
//' @docType methods
//' @description A specific z. score method that converts the data to 10 +-1 instead of 0+-1
//' in order to keep the not expressed clearly separate from the real data.
//' @param data the sparse Matrix
//' @param display_progress show a progress bar (TRUE)
//' @title Calculate z score for a sparse matrix
//' @export
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ZScore (Eigen::SparseMatrix<double> data, bool display_progress=true){


	data = data.transpose();

	Progress p(data.outerSize(), display_progress);
	//Rcout << "outerSize "<< data.outerSize() << " and innerSize " << data.innerSize() << std::endl;
	/* my normalization never creates negative values, but sores 'lost data'
	 * as -1 values. These must not be used in the normalization process
	 */
	std::vector<bool> USE(data.innerSize() );
	int iter = 0;
	//int total = 0;
	for (int k=0; k < data.outerSize(); ++k){
		//total ++;
		p.increment();
		double sum = 0.0;
		double c = 0;
		std::fill(USE.begin(), USE.end(), false );
		iter= 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if (it.value() > 0 ) {
				USE[iter] = true;
				c += 1.0;
				sum += it.value();
			}
			iter++;
		}
		double mean = sum / data.outerSize();
		//Rcout << "mean of row: "<< k << " showing a mean of " << c << " entries = "<< mean << std::endl;
		sum = 0.0;
		iter= 0;

		double sd_zero = 0.0 - mean;
		sd_zero = sd_zero * sd_zero;

		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( USE[iter++]) {
				double entry = (it.value() - mean);
				sum += entry*entry;
				it.valueRef() = entry;
				//Rcout <<  it.value() <<",";
			}
			else {
				sum += sd_zero;
			}
		}
		// sum is missing now all sd_zero parts from sparse values - add this!
		sum += sd_zero * ( data.innerSize() -c );
		double sd;


		if ( (c -1) == 0.0 ){
			sd = 0.0;
		}
		else{
			sd = sqrt(sum/( data.outerSize()-1)); //to copy the R implementation
		}

		//Rcout << " and sd = "<< sd << std::endl;
		/*Rcout << k << " mean " << mean << " and sd " << sd << "with count "<< c<< std::endl;*/
		iter= 0;
		if ( sd == 0 ){
			double entry = 10.0 ;
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
				if 	( USE[iter++]) {
					it.valueRef() = entry;
				}
			}
		}else {
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
				if 	( USE[iter++]) {
					double entry = 10.0 +(it.value() / sd)  ;
					it.valueRef() = entry;
				}
			}
		}
	}
	data = data.transpose();
	//Rcout << "Calculated a total of total different sd values: "<< total << std::endl;
	return (data);
}


//' @name ZScoreAll
//' @aliases ZScoreAll,FastWilcoxTest-method
//' @rdname ZScoreAll-methods
//' @docType methods
//' @description A sparse matrix z. score function returning the same as a normal one
//' @param data the sparse Matrix
//' @param display_progress show a progress bar (TRUE)
//' @title Calculate z score for a sparse matrix
//' @export
//[[Rcpp::export]]
NumericMatrix ZScoreAll (Eigen::SparseMatrix<double> data, bool display_progress=true){


	data = data.transpose();
	NumericMatrix res(data.outerSize(), data.innerSize());

	Progress p(data.outerSize(), display_progress);
	//Rcout << "outerSize "<< data.outerSize() << " and innerSize " << data.innerSize() << std::endl;
	/* my normalization never creates negative values, but sores 'lost data'
	 * as -1 values. These must not be used in the normalization process
	 */
	int iter = 0;
	//int total = 0;
	for (int k=0; k < data.outerSize(); ++k){
		//total ++;
		p.increment();
		double sum = 0.0;
		double c = 0;
		iter= 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			c += 1.0;
			sum += it.value();
			iter++;
		}
		double mean = sum / data.innerSize();
		Rcout << "mean of row: "<< k << " showing a mean of " << c << " entries = "<< mean << std::endl;
		sum = 0.0;
		iter= 0;

		double sd_zero = 0.0 - mean;
		sd_zero = sd_zero * sd_zero;

		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			double entry = (it.value() - mean);
			sum += entry*entry;
			it.valueRef() = entry;
		}
		// sum is missing now all sd_zero parts from sparse values - add this!
		sum += sd_zero * ( data.innerSize() -c );
		double sd;


		if ( (c -1) == 0.0 ){
			sd = 0.0;
		}
		else{
			sd = sqrt(sum/( data.innerSize()-1)); //to copy the R implementation
		}

		Rcout << " and sd = "<< sd << std::endl;
		Rcout << k << " mean " << mean << " and sd " << sd << "with count "<< c<< std::endl;
		
		iter= 0;
		
		if ( sd == 0.0 ){
			NumericVector na_vec(res.ncol(), 0.0 );
			res(k, _ ) = na_vec;
		}
		else {
			NumericVector na_vec(res.ncol(), (-mean / sd) );
			res(k, _ ) = na_vec;
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
				double entry = (it.value() / sd);
				// or vice versa ?!
				res(it.col(), it.row()) = entry;
			}
		}
	}
	data = data.transpose();
	//Rcout << "Calculated a total of total different sd values: "<< total << std::endl;
	return (res);
}

//' @name MEAN_STD
//' @aliases MEAN_STD,FastWilcoxTest-method
//' @rdname MEAN_STD-methods
//' @docType methods
//' @description A specific z. score method that converts the data to 10 +-1 instead of 0+-1
//' in order to keep the not expressed clearly separate from the real data.
//' @param data a spare matrix
//' @title Calculate mean and std of >0 values in a sparse matrix
//' @export
//[[Rcpp::export]]
NumericMatrix MEAN_STD (Eigen::SparseMatrix<double> data){
	data = data.transpose();

	NumericMatrix res(data.cols(), 4);
	std::vector<bool> USE(data.innerSize() );
	int iter = 0;
	//int total = 0;
	for (int k=0; k < data.outerSize(); ++k){
		//total ++;
		double sum = 0.0;
		int c = 0;
		std::fill(USE.begin(), USE.end(), false );
		iter= 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if (it.value() > 0 ) {
				USE[iter] = true;
				c ++;
				sum += it.value();
			}
			iter++;
		}
		double mean = sum / c;
		res(k,3) = c;
		res(k,0) = mean;
		//Rcout << "mean of row: "<< k << " showing a mean of " << c << " entries = "<< mean << std::endl;
		sum = 0.0;
		iter= 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( USE[iter++]) {
				double entry = (it.value() - mean);
				sum += entry*entry;
				it.valueRef() = entry;
				//Rcout <<  it.value() <<",";
			}
		}
		double sd = sqrt(sum/(c -1)); //to copy the R implementation

		res(k,1) = sum;
		res(k,2) = sd;
	}

	data = data.transpose();
	colnames(res) = CharacterVector::create("mean", "v-mean squared", "std", "N" );
	return res;
}

//' @name SQRT
//' @aliases SQRT,FastWilcoxTest-method
//' @rdname SQRT-methods
//' @docType methods
//' @description A specific z. score method that converts the data to 10 +-1 instead of 0+-1
//' in order to keep the not expressed clearly separate from the real data.
//' @param data a vector
//' @title Calculate z score for a sparse matrix
//' @export
//[[Rcpp::export]]
std::vector<double> SQRT (std::vector<double> data){
	std::vector<double> ret ( data.size() );
	for ( int i = 0; i < data.size(); i++ ) {
		ret[i] = sqrt( data[i] );
	}
	return ret;
}
