#' @name Rzscore
#' @aliases Rzscore,FastWilcoxTest-method
#' @rdname Rzscore-methods
#' @docType methods
#' @description calculate the z.score as the cpp function does it
#' @param obj the FastWilcoxTest
#' @title description of function Rzscore
#' @export 

setGeneric('Rzscore', ## Name
	function ( obj ) { 
		standardGeneric('Rzscore')
	}
) 

setMethod('Rzscore', signature = c ('FastWilcoxTest'),
	definition = function ( obj ) {
		i = 0
	obj@dat = Matrix::Matrix(t(
				apply(obj@dat,1, function (x) {
							i = i+1
							n <- which(x <= 0)
							dropped = which(x == -1)
							if ( length(x) - length(n) > 1 ){
								if (length(n) == 0 ){
									x <-  10 + scale(as.vector(t(x)))
								}
								else {
									x[-n] <- 10 + scale(as.vector(t(x[-n])))
									b = which(is.na(x))
									if ( length(b) > 0 ){
										x[b] = 10
									}
									#			x[n] <- 0
									#			if ( length(dropped) > 0) {
									#				x[dropped] <- -1
									#			}
								}
								
							}
							else {
								#		x[] = 0
							}
							x}
				)
		), sparse =T)
	invisible(obj)
} )


#' @name RzscoreAll
#' @aliases RzscoreAll,FastWilcoxTest-method
#' @rdname RzscoreAll-methods
#' @docType methods
#' @description calculate the z.score as the cpp function does it
#' @param obj the FastWilcoxTest
#' @title description of function Rzscore
#' @export 
setGeneric('RzscoreAll', ## Name
	function ( obj ) { 
		standardGeneric('RzscoreAll')
	}
) 

setMethod('RzscoreAll', signature = c ('FastWilcoxTest'),
	definition = function ( obj ) {
		t(apply( obj@dat, 1 , function(d) { (d - mean(d)) / sd(d) } ))
} )