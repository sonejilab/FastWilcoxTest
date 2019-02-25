#' @name print
#' @aliases print,RcppTestArea-method
#' @rdname print
#' @docType methods
#' @description  print the RcppTestArea
#' @param x the RcppTestArea object
#' @return nothing
#' @title description of function print
#' @export 
if ( ! isGeneric('print') ){setGeneric('print', ## Name
	function (x) { 
		standardGeneric('print')
	}
) }


setMethod('print', signature = c ('RcppTestArea'),
	definition = function (x) {
	cat (paste("An object of class", class(x)),"\n" )
	#cat("named ",object$name,"\n")
	cat (paste( 'with',nrow(x@dat),'genes and', ncol(x@dat),' samples.'),"\n")
	#cat (paste("Annotation datasets (",paste(dim(object$annotation),collapse=','),"): '",paste( colnames(object$annotation ), collapse="', '"),"'  ",sep='' ),"\n")
	#cat (paste("Sample annotation (",paste(dim(object$samples),collapse=','),"): '",paste( colnames(object$samples ), collapse="', '"),"'  ",sep='' ),"\n")
	#if ( length(names(object$stats)) > 0 ){
#		cat ( "P values were calculated for ", length(names(object$stats)) -1, " condition(s)\n")
	#}
})


