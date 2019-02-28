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
#' @title description of function Rstats
#' @export 
setGeneric('Rstats', ## Name
	function ( X, interest, backgound,  logFCcut = 1.0 ) { ## Argumente der generischen Funktion
		standardGeneric('Rstats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('Rstats', signature = c ('dgCMatrix'),
	definition = function ( X, interest, backgound,  logFCcut = 1.0 ) {
	logRF = unlist(apply( X, 2, function( x ) { logFC( x[interest], x[backgound] ) } ) )
	#browser()
	OK = which( logRF > logFCcut )
	W = vector( 'numeric', length(OK))
	p = vector( 'numeric', length(OK))
	id = 1;
	for ( i in OK ) {
		t = wilcox.test( X[interest, i], X[backgound,i], 'greater')
		W[id] = t$statistic
		p[id] = t$p.value
		id = id +1
	}
	
	cbind( colID = OK, logFC = logRF[OK], 'rank.sum' = W,  p.value = p)
} )

setMethod('Rstats', signature = c ('dgeMatrix'),
		definition = function ( X, interest, backgound,  logFCcut = 1.0 ) {
			logRF = unlist(apply( X, 2, function( x ) { logFC( x[interest], x[backgound] ) } ) )
			#browser()
			OK = which( logRF > logFCcut )
			W = vector( 'numeric', length(OK))
			p = vector( 'numeric', length(OK))
			id = 1;
			for ( i in OK ) {
				t = wilcox.test( X[interest, i], X[backgound,i], 'greater')
				W[id] = t$statistic
				p[id] = t$p.value
				id = id +1
			}
			
			cbind( colID = OK, logFC = logRF[OK], 'rank.sum' = W,  p.value = p)
		} )