#' @name show
#' @aliases show,FastWilcoxTest-method
#' @rdname show
#' @docType methods
#' @description  show the FastWilcoxTest
#' @param object the FastWilcoxTest object
#' @title description of function show
#' @export 
setMethod('show', signature ='FastWilcoxTest',
	definition = function (object) {
	cat (paste("An object of class", class(object)),"\n" )
	cat (paste( 'with',nrow(object@dat),'genes and', ncol(object@dat),' samples.'),"\n")
	invisible(NULL)
})


