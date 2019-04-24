/* The origin of this file was obtained from the BioQC bioconductor package (GLP-3).
 * Thanks to Jitao David Zhang, Laura Badi, Gregor Sturm and Roland Ambs, the autors of that package.
 * I pushed the original c logics into c++ classes (the header file).
 * In addition to that the main input is now either two vectors or
 * a sparse matrix instead of a conventional matrix in the original implementation.
 *
 * ! \file rank.c
  \brief statistical ranking

  Functions for statistical (fractional) ranking
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stat_rank.h"

/*! \brief Create a DRank object
 * 
 * A item object holds a double value, its original index, and its rank.
 *
 * The rank is initialized with -1, and changed to a positive one (starting from 1) by sortRankDRankList. This is used to check whether that function has been run or not, so please do not change the initial value of rank.
 */



