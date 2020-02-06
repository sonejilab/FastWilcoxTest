#include <Rcpp.h>
using namespace Rcpp;
//#include <progress.hpp> 

//' @title extract proximity for the ranger results
//' copied from https://github.com/imbs-hl/ranger/issues/234
//' @aliases extract_proximity_oob,FastWilcoxTest-method
//' @rdname extract_proximity_oob
//' @description calculate the proximity matrix
//' @param pred the predictions created from a predict(rangerRF, data, type = "terminalNodes")$predictions
//' @param prox an empty matrix with dim(pred) dimensions
//' @param inbag the inbag information from the ranger prediction run (rangerRF$inbag.counts)
//' @return prox with correct values
//' @export
// [[Rcpp::export]]
NumericMatrix extract_proximity_oob(NumericMatrix pred, NumericMatrix prox, NumericMatrix inbag) {
  
  //Progress p( prox.nrow(), true);
  // define variables
  NumericVector tree_idx;
  double same_bag_total;
  double pred_total;
  int ntree = pred.ncol();
  int ns = prox.nrow();
  
  for(int i = 0; i < ns; i++ ){
    // loop down the rows

    //p.increment();
    
    for(int j = 0; j < ns; j ++){

      if(i == j){
        // self similarity
        prox(i, j) = 1;
      }
      else{
        
        if(j < i){
          // we are below the diagonal, we can copy down the already computed similarity
          prox(i, j) = prox(j, i);
        }else{
       
        tree_idx = inbag(i, _) + inbag(j, _);
        
        pred_total = 0;
        same_bag_total = 0;
        
        for(int nt = 0; nt < ntree; nt ++){
          if( tree_idx[nt] == 0){
           
            same_bag_total += 1;
            
            if(pred(i, nt) == pred(j, nt)){
              pred_total +=  1;
            }
            
          } // end if
          
        } // end for nt
        if ( same_bag_total > 0 ){       
          prox(i, j) = pred_total / same_bag_total;
        }else{
          prox(i, j) = 0.0;
        }
        }
      }
    } // end for j
    
  } // end for i
  
  return prox;
}