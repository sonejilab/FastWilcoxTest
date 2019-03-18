
meltSparseMatrix <- function( x ){
	mdc = data.frame(
			gene_id= x@i +1, 
			cell_id = toColNums( x ), 
			value= x@x
	)
}
