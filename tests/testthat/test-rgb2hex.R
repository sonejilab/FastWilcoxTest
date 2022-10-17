context('rgb2hex')

expect_equal( c("#FFFFFF", "#FFFFFF"), rgb2hex(matrix ( rep(255, 6), nrow=2)) )

expect_equal( "#FFFFFF", rgb2hexS( 255, 255, 255 ) )

expect_equal( "#190101", rgb2hexS( 25, 2, 5) )

expect_equal( '#CEC8D2', rgb2hexS( 206, 200, 210 ) )