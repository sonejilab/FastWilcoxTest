#' @name CalcLinLang
#' @aliases CalcLinLang,FastWilcoxTest-method
#' @rdname CalcLinLang-methods
#' @docType methods
#' @description Identify genes slowly rising above the detection limit
#' @param X the sparse Matrix (row = cells, col = genes)
#' @param Grouping a numeric vector of group IDs
#' @param nGroup the number of groups
#' @param minPct ignore genes with less than 10 percent of the cells expressing them
#' @param display_progress show a progress bar (TRUE)
#' @title CalcLinLang test for rise above the detection limit
#' @export
CalcLinLang <- function(X, Grouping, nGroup, minPct = 0.1, display_progress = TRUE ){
	if ( is.factor(Grouping)) {
		## a factor needs to be converted into a vector of 1:n in the order of the factor
		new= as.vector(Grouping)
		p = names(table(Grouping)[which(table(Grouping) > 0 )])
		n = 1
		for ( i in p ) { new[which( Grouping == i ) ] = n; n=n+1}
		new = factor(new, levels= 1:length(p) )
		res = LinLang( X, new, length(p), minPct, display_progress)
	}else {
		res = LinLang( X, Grouping, nGroup, minPct, display_progress)
	}
	rownames(res) = colnames(X)
	res = cbind( res, maxR = apply( res[,1:3],1,function(x) {
						x[which(is.na(x))] = 0
						if ( x[1] * x[2] < 0 ){
							## pos and neg correlation? Not OK
							0.0
						}else {
							x[which( abs(x) == max(abs(x)))[1]] 
						}
						}))
	res
}