#' @name FastWilcoxTest
#' @title FastWilcoxTest
#' @description  An S4 class to train to code R C++
#' The class itself is totally unused. All function are directly applied to the dgCMatrix matrix.
#' @slot dat a dgCMatrix to learn how to modify that using C++
#' @slot usedObj here a set of used and probably important objects can be saved. Be very careful using any of them!
#' @exportClass FastWilcoxTest
setClass(
		Class='FastWilcoxTest', 
		representation = representation ( 
				dat='dgCMatrix',
				usedObj= 'list'
		),
		prototype(
				dat=Matrix::Matrix(rep(0,2),sparse=T),
				usedObj= list()
  )
)