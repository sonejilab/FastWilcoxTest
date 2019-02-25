
Rzscore <- function ( obj ) {
	obj@dat = Matrix(t(
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
		))
	invisible(obj)
}
