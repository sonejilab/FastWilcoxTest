#' @name Rstats
#' @aliases Rstats,FastWilcoxTest-method
#' @rdname Rstats-methods
#' @docType methods
#' @description calculate the wilcox test on a sparse matrix as the C++ function does it.
#' The data in the sparse matrix should be log transformed (not log10 or log2!)
#' @param X the sparse matrix
#' @param interest the intereting col IDs
#' @param backgound the background col IDs
#' @param logFCcut the logFC change cutoff default= 1.0
#' @param minPct only test genes that are detected in a minimum fraction of
#'  min.pct cells in either of the two populations. Meant to speed up the function
#' @param onlyPos test only those genes with higher expression in the interest (default FALSE)
#' @title description of function Rstats
#' @export 
setGeneric('Rstats', ## Name
	function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1 , onlyPos=FALSE) { ## Argumente der generischen Funktion
		standardGeneric('Rstats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('Rstats', signature = c ('dgCMatrix'),
	definition = function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1, onlyPos=FALSE ) {
		Rstats( as.matrix(X), interest, backgound,  logFCcut, minPct, onlyPos )
} )

setMethod('Rstats', signature = c ('matrix'),
		definition = function ( X, interest, backgound,  logFCcut = 1.0, minPct=0.1, onlyPos=FALSE ) {
			logRF = unlist(apply( X, 2, function( x ) { logFC( x[interest], x[backgound] ) } ) )
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
				t = wilcox.test( X[interest, i], X[backgound,i])
				W[id] = t$statistic
				p[id] = t$p.value
				id = id +1
			}
			
			cbind( colID = OK, logFC = logRF[OK], 'rank.sum' = W,  p.value = p)
		} )