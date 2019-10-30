#' @name euclidian_order3d
#' @aliases euclidian_order3d,FastWilcoxTest-method
#' @rdname euclidian_order3d-methods
#' @docType methods
#' @description identify a order in a linear set of 3D coordinates.
#' @param x  the x coordinate
#' @param y  the y coordinate
#' @param z  the z coordinate
#' @param last not re-sample the data (default = FALSE)
#' @title description of function euclidian_order3d
#' @export 
setGeneric('euclidian_order3d', ## Name
	function ( x, y, z, last=FALSE ) { ## Argumente der generischen Funktion
		standardGeneric('euclidian_order3d') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('euclidian_order3d', signature = c ('numeric'),
	definition = function ( x, y, z, last=FALSE ) {

	o = eDist3d( x, y, z ,0)
	order = vector('numeric', length(x))
	order[1] = order(o)[length(o)]
	origIDs = 1:length(o)
	for( i in 2:length(x)) { 
		OK = origIDs
		a = i-1
		if ( i >3 ) {
			OK = origIDs[ - order[1:(i-2)] ]
		}
		if ( length(OK) == 1) {
			order[i] = origIDs[OK]
			}else {
			if ( is.na( match( order[a], OK ) ) ) {
				## OK this was a bad starting cell - re-sample that
				change = sample( 1:length(x)) 
				ot =  euclidian_order3d( x[change], y[change], z[change] , TRUE)
				return( change[ot] )
			}
			t =order( eDist3d( x[OK], y[OK], z[OK] , match( order[a], OK ) -1 ) )
			order[i] = OK[t[2]]
		}
	}
	order
} )
