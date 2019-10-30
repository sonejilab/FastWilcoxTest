
euclidian_order3d <- function( x, y, z ){

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
				browser()
			}
			t =order( eDist3d( x[OK], y[OK], z[OK] , match( order[a], OK ) -1 ) )
			order[i] = OK[t[2]]
		}
	}
	order
}
