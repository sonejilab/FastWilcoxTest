
if ( ! isGeneric('as_FastWilcoxTest') ){setGeneric('as_FastWilcoxTest', ## Name
	function ( x ) { 
		standardGeneric('as_FastWilcoxTest')
	}
) }

#' @name as_FastWilcoxTest
#' @aliases as_FastWilcoxTest,matrix-method
#' @rdname as_FastWilcoxTest
#' @docType methods
#' @description Creates a FastWilcoxTest from a matrix.
#' @param x a matrix object
#' @title description of function as_FastWilcoxTest
#' @examples
#' set.seed(1)
#' ncol = 100
#' nrow=90
#' dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
#' dat[which(dat < 1)] = 0
#' colnames(dat) <- paste('Sample', 1:ncol)
#' rownames(dat) <- paste( 'gene', 1:nrow)
#' 
#' x <- as_FastWilcoxTest( dat )
#' x
#' @export 
setMethod('as_FastWilcoxTest', signature = c ('matrix'),
	definition = function ( x ) {
	methods::new('FastWilcoxTest', dat=Matrix::Matrix( x, sparse=T) )
} )
