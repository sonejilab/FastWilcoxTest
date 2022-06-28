setGeneric('Rstats', ## Name
	function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1 , onlyPos=FALSE) { ## Argumente der generischen Funktion
		standardGeneric('Rstats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

#' @name Rstats
#' @aliases Rstats,dgCMatrix-method
#' @rdname Rstats-methods
#' @docType methods
#' @description calculate the wilcox test on a sparse matrix as the C++ function does it.
#' This function has been implemented for test puroses and keeped to allow tests for every user.
#' @param X the sparse matrix
#' @param interest the interesting col IDs
#' @param backgound the background col IDs
#' @param logFCcut the logFC change cutoff default= 1.0
#' @param minPct only test genes that are detected in a minimum fraction of
#'  min.pct cells in either of the two populations. Meant to speed up the function
#' @param onlyPos test only those genes with higher expression in the interest (default FALSE)
#' @title Run the 'wilcoxTest' on a sparse matrix like the c++ 'StatsTest'
#' @export 
setMethod('Rstats', signature = c ('dgCMatrix'),
	definition = function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1, onlyPos=FALSE ) {
		Rstats( as.matrix(X), interest, backgound,  logFCcut, minPct, onlyPos )
} )

#' @describeIn Rstats matrix
#' @docType methods
#' @description use the R wilcox.test function instead of a c++ version.
#' @param X the normal matrix
#' @param interest the interesting col IDs
#' @param backgound the background col IDs
#' @param logFCcut the logFC change cutoff default= 1.0
#' @param minPct only test genes that are detected in a minimum fraction of
#'  min.pct cells in either of the two populations. Meant to speed up the function
#' @param onlyPos test only those genes with higher expression in the interest (default FALSE)
#' @title description of function Rstats
#' @export
setMethod('Rstats', signature = c ('matrix'),
		definition = function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1, onlyPos=FALSE ) {
			logRF = unlist(apply( X, 2, function( x ) { FC( x[interest], x[backgound] ) } ) )
			fracA = unlist(apply( X, 2, function( x ) { x =  x[interest]; length( x[which(x > 0)] ) /length(x) } ) )
			fracB = unlist(apply( X, 2, function( x ) { x =  x[backgound]; length( x[which(x > 0)] ) /length(x) } ) )
			#browser()
			OK = NULL
			if ( onlyPos ) {
				OK = which( logRF > logFCcut )
			}else {
				OK = which( abs(logRF) > logFCcut )
			}
			tmp = unique( c( which( fracA[OK] > minPct)),names(which( fracB[OK] > minPct) ) )

			if ( length(tmp) == 0 ) {
				stop( "No genes passed the cutoff" )
			}
			OK = OK[tmp]
			W = vector( 'numeric', length(OK))
			p = vector( 'numeric', length(OK))
			id = 1;
			message ( paste( length(OK), "genes pass the logFC and frac. expr. filters."))
			for ( i in OK ) {
				t = stats::wilcox.test( X[interest, i], X[backgound,i], exact = FALSE)
				W[id] = t$statistic
				p[id] = t$p.value
				id = id +1
			}
			
			cbind( colID = OK, logFC = logRF[OK], 'rank.sum' = W,  p.value = p)
		} )
