#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

#ifndef __APPLE__
#include "omp.h"
#endif

#include "stat_rank.h"
/*#include "wmw_test.h"*/

#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#define ABSLOG(x) fabs(log10( (x) ))

// [[Rcpp::interfaces(r, cpp)]]

/*#include <Rcpp.h>

using namespace Rcpp;*/

/*	greater=0,
	less=1,
	twoSided=2,
	U=3,
	abslog10greater=4,
	log10less=5,
	abslog10twoSided=6,
	Q=7
*/

/*
 * void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)

 * i_tail in {0,1,2} means: "lower", "upper", or "both" :
 * if(lower) return  *cum := P[X <= x]
 * if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
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
      pnorm_both(zval, &pgt, &plt, 0, 0);
      res = type==greater ? pgt : ABSLOG(pgt);
    } else if (type == 1 || type == 5) { /* less */
zval = (uStat-0.5-mu)/sqrt(sigma2); // z higher tail
      pnorm_both(zval, &pgt, &plt, 1, 0);
      res = type==less ? plt : log10(plt);
    } else if (type == 2 || type == 6 || type == 7) { /* two sided*/
zval = (uStat-mu-(uStat>mu ? 0.5 : -0.5))/sqrt(sigma2);
      pnorm_both(zval, &pgt, &plt, 2, 0);
      res = mu==0.0 ? 1.0 : 2.0*MIN(pgt, plt);
      if(type == 4) {
        res = ABSLOG(res);
      } else if (type == 7) {
        res = pgt<=plt ? ABSLOG(res) : -ABSLOG(res);
      }
    } else {
      error("Unrecognized type %d. Should not happen\nPossible only int values  0=greater, 1=less, 2=twoSided, 3=U, 4=abslog10greater, 5=log10less, 6=abslog10twoSided, 7=Q",
            type);
    }
  }
  return(res);
}
/*
 double wmw_test_core(const DRankList valList,
 const int *inds, int nInds,
 int nTotal, TestType type) {
 int i;
 double indRankSum; // sum of index rank
 double res;
 
 indRankSum = 0.0;
 for(i = 0;i<nInds;++i)
 indRankSum += valList->list[inds[i]]->rank;
 
 res = wmw_test_stat(indRankSum, nInds, nTotal,
 tieCoef(valList), type);
 return(res);
 }
 */

void wmw_test_list(const double *valPtr, int n,
                   SEXP indlist,
                   double *resPtr, int type) {
  DRankList list;
  int i, j;
  int nInd;
  int* ip;
  
  double tie;
  double indRankSum;
  
  list = createDRankList(valPtr, n);
  prepareDRankList(list);
  
  tie = tieCoef(list);
  
#pragma omp parallel for
  for(i=0;i<length(indlist);++i) {
    ip=INTEGER(VECTOR_ELT(indlist,i));
    nInd=length(VECTOR_ELT(indlist,i));
    
    indRankSum = 0.0;
    for(j=0; j<nInd; ++j) {
      if(!(ip[j]>=0 && ip[j]<=n-1))
        error("Index out of range: gene set %d, gene %d\n", i+1, j+1);
      indRankSum += list->list[ip[j]]->rank;
    }
    
    resPtr[i] = wmw_test_stat(indRankSum,
                              nInd,
                              n,
                              tie,
                              type);
  }
  destroyDRankList(list);
}


/*! \brief Wilcoxon-Mann-Whitney Test
 *
 * \param indlist: A list of integers giving the index of gene sets
 * \param matrix: an expression matrix with features in rows and samples in columns
 * \param rtype:
 * \parblock
 * Define f(x)=abs(log10(x))
 * 0=p(greater), 1=p(less), 2=p(twoSided), 3=U,
 * 4=f(p(greater)),5=p(less),6=f(p(twoSided)), 7=p(greater)<p(less) ? f(p(twoSided)) : -f(p(twoSided))
 * \endparblock
 *
 * This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
 */
extern SEXP wmw_test(const double vec, SEXP indlist, int rtype) {
/* const int type=INTEGER(rtype)[0]; */
  const int m=length(indlist);
  const int n=length(vec);
  
  int i;
  double *matColPtr; // pointer to the current column of the matrix
  SEXP res;
  double *resPtr;
  
  res=PROTECT(allocMatrix(REALSXP, m, NCOL(matrix)));
  
  resPtr=REAL(res);
  
  wmw_test_list(vec, n, indlist, resPtr, rtype);
  
  return(res);
}

// ----------------------------------------
// signed tests
// ----------------------------------------
void signed_wmw_test_list(const double *valPtr, int n,
                          SEXP signedIndList,
                          double *resPtr, int type) {
  DRankList list;
  int i, j;
  int nPos, nNeg;
  int *ipPos,*ipNeg;
  SEXP pairInd, posInd, negInd;
  
  double tie;
  double indRankSum;
  
  list=createDRankList(valPtr, n);
  prepareDRankList(list);
  
  tie = tieCoef(list);
  
#pragma omp parallel for
  for(i=0;i<length(signedIndList);++i) {
    pairInd = VECTOR_ELT(signedIndList, i);
    posInd = VECTOR_ELT(pairInd, 0);
    negInd = VECTOR_ELT(pairInd, 1);
    
    indRankSum = 0.0;
    
    if(posInd != NULL_USER_OBJECT) {
      ipPos = INTEGER(posInd);
      nPos = length(posInd);
      for(j=0; j<nPos; ++j) {
        if(!(ipPos[j]>=0 && ipPos[j]<=n-1))
          error("Index out of range: gene set %d, gene %d\n", i+1, j+1);
        indRankSum += list->list[ipPos[j]]->rank;
      }
    } else {
      nPos = 0;
    }
    
    if(negInd != NULL_USER_OBJECT) {
      ipNeg = INTEGER(negInd);
      nNeg = length(negInd);
      for(j=0; j<nNeg; ++j) {
        if(!(ipNeg[j]>=0 && ipNeg[j]<=n-1))
          error("Index out of range: gene set %d, gene %d\n", i+1, j+1);
        indRankSum += (n - list->list[ipNeg[j]]->rank + 1);
      }
    } else {
      nNeg = 0;
    }
    
    resPtr[i]=wmw_test_stat(indRankSum,
                            nPos+nNeg,
                            n,
                            tie,
                            type);
  }
  destroyDRankList(list);
}

/*! \brief Signed Wilcoxon-Mann-Whitney Test
 *
 * \param matrix: an expression matrix with features in rows and samples in columns
 * \param signedIndList: A list of signed gene-set-pairs, each of which is a list of two character vectors, containing positive and negative signatures, respectively
 * \param rtype:
 * \parblock
 * Define f(x)=abs(log10(x))
 * 0=p(greater), 1=p(less), 2=p(twoSided), 3=U,
 * 4=f(p(greater)),5=p(less),6=f(p(twoSided)), 7=p(greater)<p(less) ? f(p(twoSided)) : -f(p(twoSided))
 * \endparblock
 *
 * This implementation uses normal approximation, which works reasonably well if sample size is large (say N>=20)
 */

extern SEXP  signed_wmw_test(SEXP matrix, SEXP signedIndList, int rtype) {
  /*const int type=INTEGER(rtype)[0];*/
  const int m=length(signedIndList);
  const int n=NROW(matrix);
  
  int i;
  double *matColPtr; // pointer to the current column of the matrix
  SEXP res;
  double *resPtr;
  
  res=PROTECT(allocMatrix(REALSXP, m, NCOL(matrix)));
  
  resPtr=REAL(res);
  matColPtr=REAL(matrix);
  
#pragma omp parallel for
  for(i=0; i<NCOL(matrix);++i) {
/*cout << i << " n " << n << " signedIndList " << signedIndList << " resPtr "<< resPtr << std::endl; */

	  signed_wmw_test_list(matColPtr, n,
                         signedIndList,
                         resPtr, rtype);
    resPtr+=m;
    matColPtr+=n;
  }
  
  UNPROTECT(1);
  return( res );
}
