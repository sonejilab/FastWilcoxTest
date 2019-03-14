// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
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
	/* Progress p(data.outerSize(), display_progress); */
	data = data.transpose();

	//int total = 0;
	for (int k=0; k < data.outerSize(); ++k){
		//total ++;
		/*p.increment();*/
		double sum = 0.0;
		double c = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			c += 1.0;
			sum += it.value();
		}
		double mean = sum / c;
		//Rcout << "mean of row: "<< k << " showing a mean of " << c << " entries = "<< mean << std::endl;
		sum = 0.0d;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			double entry = (it.value() - mean);
			sum += entry*entry;
			it.valueRef() = entry;
			//Rcout <<  it.value() <<",";
		}
		double sd = sqrt(sum/c);

		//Rcout << " and sd = "<< sd << std::endl;
		/*Rcout << k << " mean " << mean << " and sd " << sd << "with count "<< c<< std::endl;*/
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			double entry = 10.0 +(it.value() / sd)  ;
			it.valueRef() = entry;
		}
	}
	data = data.transpose();
	//Rcout << "Calculated a total of total different sd values: "<< total << std::endl;
	return (data);
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

	NumericMatrix res(data.cols(), 3);
	for (int k=0; k < data.outerSize(); ++k){
		double sum = 0.0;
		double c = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( it.value() != -1) {
				c += 1.0;
				sum += it.value();
			}
		}
		double mean = sum / c;
		res(k,0) = mean;
		//Rcout << "mean of row: "<< k << " showing a mean of " << c << " entries = "<< mean << std::endl;
		sum = 0.0d;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			if ( it.value() != -1) {
				double entry = (it.value() - mean);
				sum += entry*entry;
				it.valueRef() = entry;
			}
			//Rcout <<  it.value() <<",";
		}
		double sd = sqrt(sum/c);
		res(k,1) = sum;
		res(k,2) = sd;
	}

	data = data.transpose();
	colnames(res) = CharacterVector::create("mean", "v-mean squared", "std" );
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
