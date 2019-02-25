#' @name as_RcppTestArea
#' @aliases as_RcppTestArea,RcppTestArea-method
#' @rdname as_RcppTestArea-methods
#' @docType methods
#' @description Creates a RcppTestArea from a matrix.
#' @param x a matrix object
#' @title description of function as_RcppTestArea
#' @export 
if ( ! isGeneric('as_RcppTestArea') ){setGeneric('as_RcppTestArea', ## Name
	function ( x ) { 
		standardGeneric('as_RcppTestArea')
	}
) }

setMethod('as_RcppTestArea', signature = c ('matrix'),
	definition = function ( x ) {
	new('RcppTestArea', dat=Matrix::Matrix( x, sparse=T) )
} )
