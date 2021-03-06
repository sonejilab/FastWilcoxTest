context( 'Cpp euclidian order')

prefix='./'
# prefix='tests/testthat'; devtools::load_all()

data = readRDS( file.path( prefix, 'data', 'problematixTimeline.rds'))


o = eDist3d( data$x, data$y, data$z ,0)
order = vector('numeric', length(data$x))
order[1] = order(o)[length(o)]
origIDs = 1:length(o)


for( i in 2:length(data$x)) { 
		OK = origIDs
		a = i-1
		if ( i >2 ) {
			OK = origIDs[ - order[1:(i-2)] ]
		}
		if ( length(OK) == 1) {
			order[i] = origIDs[OK]
		}
		else {
			if ( is.na( match( order[a], OK ) ) ) {
				message('OK this was a bad starting cell - re-sample that')
				change = sample( 1:length(data$x)) 
				if ( last ) {
					browser()
				}else {
					ot =  euclidian_order3d( data$x[change], data$y[change], data$z[change] , TRUE)
				}
				return( change[ot] )
			}
			t =order( eDist3d( data$x[OK], data$y[OK], data$z[OK] , match( order[a], OK ) -1 ) )
			order[i] = OK[t[2]]
		}
	}

expect_equal( euclidian_order3d( data$x, data$y, data$z, FALSE ), order )


for ( i in 1:300 ) {
	change = sample( 1:length(data$x)) 
	cat('.')
	ot =  euclidian_order3d( data$x[change], data$y[change], data$z[change], FALSE )
	if ( !  all.equal( change[ot], order) == TRUE ) {
		ot = rev(ot)
	}
	expect_equal(change[ot], order)
}

#rgl::plot3d( data$x[order], data$y[order], data$z[order], col=gplots::bluered( length(data$x)))

#o = euclidian_order3d( data$x, data$y, data$z )