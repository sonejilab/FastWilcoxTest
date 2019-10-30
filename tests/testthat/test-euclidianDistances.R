context( 'Cpp euclidian disnatces')

prefix='./'

x=  c( 1:10 )
y = c ( 1:10 )

dist = sqrt( 2 )

result = euclidian_distances( x, y )
exp = c( 0, rep( dist, 9))

expect_equal( result, exp, 1e-7)

