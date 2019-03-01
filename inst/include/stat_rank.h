/*! \file stat_rank.h
    \brief Header for statistical ranking

    Details.
 */

#include <Rcpp.h>
using namespace Rcpp;



class DRank {
public:
	int index; /*!< Input index (starting from 0) */
	double vPtr; /*!< Pointer to the value (to use with R) */
	double rank; /*!< Fractional ranking (starting from 1) */
	DRank () { index = -1; vPtr = -1; rank = -1; };
	DRank (int ind, double Ptr) { index = ind; vPtr = Ptr; rank = -1; };
	void fill (int ind, double Ptr ) {
		index = ind;
		vPtr = Ptr;
		rank =-1;
	};
	void setRank( double r ) { if ( r > 0 ) { rank = r; } };
	void print( void ){
		Rcout << "C++ DRank entry at index " << index << ", rank " << rank << " with value "<< vPtr << std::endl;
	};
};

/* following this:
 * https://stackoverflow.com/questions/20158793/creating-c-vector-of-pointers
 * https://de.cppreference.com/w/cpp/memory/unique_ptr
 * http://www.cplusplus.com/reference/vector/vector/push_back/
 * https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique
 * */

class DRankList {
public:

	std::vector<DRank> list; /*!< Dynamic array of DRanks */
	int len; /*!< Length of the array */
	int ulen; /*!< Length of unique elements */
	double tieCoef; /*!< Tie coefficient used by the WMW test */

	/*! \brief create an DRankList object
     \param array: a double array
     \param len: the length of the double array
	 */
	DRankList() { len = 0; tieCoef = -1; ulen = 0;};
	void refill( std::vector<double> array, int newlen) {
		// in case we populate first time
		//Rcout << "C++ DRankList you ask for a refill ;-) " << newlen << " == my length? " << len << std::endl;
		if ( len == 0 ) {
			len = newlen;
			//Rcout << "C++ DRankList you ask for a refill ;-) I try to create the vector" << std::endl;
			list.clear(); // to be save here!
			for(int i=0;i<len;++i) {
				list.push_back( DRank (  i, array.at(i)) );
			}
			//Rcout << "and it worked!" << list.size() << " is the new size" << std::endl;
			return;
		}else if ( len ==  newlen) {
			// just re_init them later
			//Rcout << "using the old data as length should be ok" << std::endl;
		}else { // new length!!
			//Rcout << "resetting the vector to new length" << std::endl;
			len = newlen;
			list.clear();
			list.resize(len);
		}
		ulen=-1;
		tieCoef = -1;
		for(int i=0;i<len;++i) {
			//Rcout << "C++ DRankList at pos " << i << " fill with " << array.at(i) << "and list size" << list.size() << std::endl;
			if ( i >= array.size() ) {
				::Rf_error("DRankList is bigger than data vector");
			}
			list.at(i).fill( i, array.at(i) );
		}
	};
	void print(void){
		Rcout << "C++ DRankList with " << len << " entries, " << ulen << " unique entries and "<< tieCoef << " ties:" << std::endl;
		for(int i=0;i<len;++i) {
			list.at(i).print();
		}
	};
	/*! \brief test whether the DRankList has been ranked
	 *
	 * if sortRankDRankList has been run, the value will be 1, otherwise 0.
	 */
	int isRanked(void) {
		if ( list.size() == 0)
			return 0;
		return(list.at(0).rank>0);
	};
	/*! \brief: sort and gives rank to a DRankList
	 *
	 *  sortRankDRankList sorts and gives statistical (fractional) rank to a DRankList.
	 *
	 *  sortRankDRankList runs once and only once (controlled by isRanked):
	 *  Once the ranks have been set (i.e. ranks>0), the function will exit
	 *  without doing anything.
	 */
	void sortRankDRankList() {
		if(this->isRanked()) return; // make sure that this function runs only once
		//Rcout << "C++ DRankList sortRankDRankList CALCULATING ONCE" << std::endl;
		std::vector<double> backup(len);

		int i, j, k;
		int ucount=0;

		for(i=0;i<len; ++i)
			backup.at(i)=list.at(i).vPtr;

		std::sort(list.begin(), list.end(),  [](DRank const &l, DRank const &r) { return l.vPtr < r.vPtr; } );

		for(i=0; i<len;i=j+1) {
		    j=i;
		    while((j<len-1) && backup.at(list.at(j).index) == backup.at(list.at(j+1).index) ) {
		      j++;
		    }
		    for(k=i;k<=j;k++)
		      list.at(k).rank=(i+j+2)*0.5;
		    ucount++;
		 }

		ulen=ucount;
	};
	/*! \brief: rankDRankList
	 * \param list A DRankList object
	 * It calls sortRankDRankList if the DRankList has not been ranked before
	 * The items in the list are sorted by input index
	 */
	void rankDRankList() {
		this->sortRankDRankList();
		std::sort(list.begin(), list.end(),  [](DRank const &l, DRank const &r) { return l.index < r.index; } );
	};
	/*! \brief: sortDRankList
	 * \param list A DRankList object
	 * It calls sortRankDRankList if the DRankList has not been ranked before
	 * The items in the list are sorted by ascending order of the values.
	 */
	void sortDRankList() {
		this->sortRankDRankList();
		std::sort(list.begin(), list.end(),  [](DRank const &l, DRank const &r) { return l.vPtr < r.vPtr; } );
	};

	/*! \brief: prepareDRankList
	 * \param list A DRankList object
	 * It prepares a DRankList object to be used in Wilcoxon-Mann-Whitney tests
	 */
	void prepareDRankList() {
		this->rankDRankList();
		//Rcout << "rankDRankList(ed) " << std::endl;
		//print();
		if( len == ulen ) {
			//Rcout << "C++ DRankList prepareDRankList setting tieCoef to 1.0 " << std::endl;
			//print();
			tieCoef=1.0;
		} else {
			int n=len;
			int un=ulen;
			int ncount=0;
			double mTieCoef=0.0;

			std::vector<int> tbl(ulen);

			int i, j;
			this->sortDRankList();
			for(i=0;i<n;i=j+1) {
				j=i;
				while(j<n-1 && ( list.at(j+1).vPtr ==list.at(j).vPtr ) ) ++j;
				tbl.at(ncount++)=j-i+1;
			}
			for(i=0;i<un;++i)
				mTieCoef+=(0.0+tbl.at(i))/n*(tbl.at(i)+1)/(n+1)*(tbl.at(i)-1)/(n-1);
			tieCoef=1-mTieCoef;
			this->rankDRankList();
		}
		//Rcout << "returned(ed) " << std::endl;
		//print();
	};


};

