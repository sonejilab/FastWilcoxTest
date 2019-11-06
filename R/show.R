#' @name show
#' @aliases show,BioData-method
#' @rdname show-methods
#' @docType methods
#' @description  show the BioData
#' @param x the BioData object
#' @param object  TEXT MISSING
#' @return nothing
#' @title description of function show
#' @export 
setGeneric('show', ## Name
	function (object) { ## Argumente der generischen Funktion
		standardGeneric('show') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('show', signature = c ('FastWilcoxTest'),
	definition = function (object) {
	cat (paste("An object of class", class(object)),"\n" )

	cat (paste( 'with',nrow(object@dat),'genes and', ncol(object@dat),' samples.'),"\n")
	invisible(NULL)
})

