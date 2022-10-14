context('rgb2hex')

expect_equal( c("#ffffff", "#ffffff"), rgb2hex(matrix ( rep(255, 6), nrow=2)) )

expect_equal( "ffffff", rgb2hexS( 255, 255, 255, FALSE) )

expect_equal( "190205", rgb2hexS( 25, 2, 5, FALSE) )