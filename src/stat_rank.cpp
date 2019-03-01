/*! \file rank.c
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



