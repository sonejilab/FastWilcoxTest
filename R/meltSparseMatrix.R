#' @name meltSparseMatrix
#' @aliases meltSparseMatrix,FastWilcoxTest-method
#' @rdname meltSparseMatrix-methods
#' @docType methods
#' @description melt a sparse matrix into row ids, col ids and value
#' @param x a sparse matrix
#' @title description of function meltSparseMatrix
#' @return a data.frame with the columns gene_id, cell_id and value
#' containing the row ids, the col ids and the value stored in the sparse matrix using c++
#' @export 
if ( ! isGeneric('meltSparseMatrix') ){setGeneric('meltSparseMatrix', ## Name
	function ( x ) { 
		standardGeneric('meltSparseMatrix')
	}
) }

setMethod('meltSparseMatrix', signature = c ('dgCMatrix'),
	definition = function ( x ) {
	data.frame(
			gene_id= x@i +1, 
			cell_id = toColNums( x ), 
			value= x@x
	)
} )

