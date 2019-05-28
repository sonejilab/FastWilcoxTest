

context( 'fast cor')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

t1 = rnorm( 30 )

cor = correlationCoefficient( X= t1, Y = t1 )

expect_equal( cor ,1 )

cor = correlationCoefficient( X= t1, Y =rev( t1 ))
expect_equal( cor , -0.05584, tolerance=1e-5 )

system.time({ cor = CorMatrix (x@dat, x@dat[1,] ) })

system.time({cor2 = apply(as.matrix(x@dat),1, cor, x@dat[1,] ) })


all.equal(cor, as.vector(cor2), 1e-7)


