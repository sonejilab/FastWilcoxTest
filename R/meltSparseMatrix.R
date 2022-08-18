
setGeneric('meltSparseMatrix', ## Name
	function ( x ) { 
		standardGeneric('meltSparseMatrix')
	}
) 



#' @name meltSparseMatrix
#' @aliases meltSparseMatrix,dgCMatrix-method
#' @rdname meltSparseMatrix-methods
#' @docType methods
#' @description melt a sparse matrix into a data.frame with gene.id = row ids, cell_id = col ids and value = x@x
#' @param x a sparse matrix
#' @title fast melt a sparse matrix omitting all '0' entries.
#' @export 
setMethod('meltSparseMatrix', signature = c ('dgCMatrix'),
	definition = function ( x ) {
	data.frame(
			gene_id= x@i +1, 
			cell_id = toColNums( x ), 
			value= x@x
	)
} )


# #' Popultes an SQLite table with the data. 
# #' Likely slower than to melt the data and use the RSQLite functionality instead.
# #' But more memory efficient.
# #' 
# #' @name sparse2SQLite_text_file
# #' @aliases sparse2SQLite_text_file,dgCMatrix-method
# #' @rdname sparse2SQLite_text_file-methods
# #' @docType methods
# #' @description melt a sparse matrix into a data.frame with gene.id = row ids, cell_id = col ids and value = x@x
# #' @param x a sparse matrix
# #' @title fast melt a sparse matrix omitting all '0' entries.
# #' @export 
# setMethod('sparse2SQLite_text_file', signature = c ('dgCMatrix'),
# 	definition = function ( x ) {

# 	data.frame(
# 			gene_id= x@i +1, 
# 			cell_id = toColNums( x ), 
# 			value= x@x
# 	)
# } )