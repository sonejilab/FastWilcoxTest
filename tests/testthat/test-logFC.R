
context( 'logFC')


expect_equal( logFC( 2, 8 ), log(4), label= "simple logFC 2->8" )
expect_equal( logFC( 8, 2 ), log(1/4), label= "simple logFC 8->2" )


expect_equal( logFC( c(8,8,8,8), c(2,2,2) ), log(1/4), label= "simple logFC 8->2" )

