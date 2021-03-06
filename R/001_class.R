#' The FastWilcoxTest class
#' 
#' is only a wrapper around a sparse matrix. The class itself is not used apart from the test scripts.
#' All functionallity is exported as c++ header FastWilcoxTest.h 
#' as well as normal R function FastWilcosTest::<function name>().
#' 
#' @name FastWilcoxTest-class
#' @rdname FastWilcoxTest-class
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

#' @name x
#' @title A simple FastWilcoxTest object containing random data to run the tests on.
#' @description The data was created following the example in as_FastWilcoxTest()
#' @docType data
#' @usage x
#' @format FastWilcoxTest
#' @keywords data
'x'
NULL;

#todo
