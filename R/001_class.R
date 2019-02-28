#' @name FastWilcoxTest
#' @title FastWilcoxTest
#' @description  An S4 class to train to code R C++
#' @slot dat a dgCMatrix to lear how to modify that using C++
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
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