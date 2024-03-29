
setGeneric('as_FastWilcoxTest', ## Name
	function ( x ) { 
		standardGeneric('as_FastWilcoxTest')
	}
) 


#' @name as_FastWilcoxTest
#' @aliases as_FastWilcoxTest,matrix-method
#' @rdname as_FastWilcoxTest
#' @docType methods
#' @description Creates a FastWilcoxTest from a matrix.
#' @param x a matrix object
#' @title description of function as_FastWilcoxTest
#' @export 
setMethod('as_FastWilcoxTest', signature = c ('matrix'),
	definition = function ( x ) {
	methods::new('FastWilcoxTest', dat=Matrix::Matrix( x, sparse=T) )
} )
