

#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>       /* log */
#include <vector>

class Table {
public:
	std::vector<int> size;
	std::map<double, int> name;

	Table () {  } ;
	Table( std::vector<double> data ) {
		// if this object has been used before null it
		if (size.size() > 0 ) {
			reset();
		}
		for ( unsigned i = 0; i < data.size(); i ++ ){
			add(data.at(i));
		}
	};
	void reset() {
		std::fill(size.begin(), size.end(), 0 );
	}
	void add( double d ) {
		if ( this->name.count( d ) == 0 ){
			int s = this->size.size();
			this->name.insert( std::pair<double,int>( d, s ) );
			this->size.push_back( 0 );
			 // Rcout << "creating key value pair: " <<  val << ", " << s << std::endl;
		}
		//Rcout << "adding key " <<  data.at(i) << " on position " << this->name.find( data.at(i)) ->second  << std::endl;
		this->size[this->name.find( d ) ->second ] ++;
	}
	unsigned sum() {
		unsigned sum = 0;
		for(auto& x : this->name){
			sum += this->size.at(x.second);
		}
		return (sum);
	}
	double _entropy(){
		double ret = 0.0;

		if ( this->size.size() == 0 ) {
			Rcout << "Please initialize the Table first using either Table( std::vector<doubles> ) or Entropy( std::vector<double> )" << std::endl;
		}else {
			int sum = 0;
			int not0 = 0;
			double p;

			for(auto& x : this->name){
				if ( this->size.at(x.second) != 0) {
					not0++;
					sum += this->size.at(x.second);
				}
			}
			//Rcout << "total not0 " <<  not0 << " with total sum " << sum << std::endl;
			if ( not0 > 1 ) {
				for(auto& x : this->name){
					if ( this->size.at(x.second) != 0){
						p = (double)this->size.at(x.second) / (double)sum ;
						//Rcout << "key " <<  x.first << " with value " << this->size.at(x.second) << " became p" << p << " and log(p) = "<< log( p ) << std::endl;
						ret += p * log( p );
					}	
				}
			}
		}
		//Rcout << "total enropy " << -ret << std::endl;
		return ( -ret );
	};
	double Entropy( std::vector<double> data ){
		Table Table( data );
		//Table.print();
		return (Table._entropy() );
	};
	void print( void ){
		
		for(auto& x : this->name){
			Rcout << x.first << " => " << this->size.at(x.second) << std::endl;
		}

		//Rcout << "C++ Table for Entropy calculation with" << this->size.size() << " different integer slots " << std::endl;
	};
};