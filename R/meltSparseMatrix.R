if ( ! isGeneric('meltSparseMatrix') ){setGeneric('meltSparseMatrix', ## Name
	function ( x ) { 
		standardGeneric('meltSparseMatrix')
	}
) }
#' @name meltSparseMatrix
#' @aliases meltSparseMatrix,dgCMatrix-method
#' @rdname meltSparseMatrix-methods
#' @docType methods
#' @description melt a sparse matrix into a data.frame with gene.id = row ids, cell_id = col ids and value = x@x
#' @param x a sparse matrix
#' @title fast melt a sparse matrix omitting all '0' entries.
#' @examples
#' x
#' melted = meltSparseMatrix(x@dat)
#' dim(melted)
#' colnames(melted)
#' @export 
setMethod('meltSparseMatrix', signature = c ('dgCMatrix'),
	definition = function ( x ) {
	data.frame(
			gene_id= x@i +1, 
			cell_id = toColNums( x ), 
			value= x@x
	)
} )

