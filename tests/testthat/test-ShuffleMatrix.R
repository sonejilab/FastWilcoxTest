set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

d = ShuffleMatrix( x@dat)

t1 =apply( d, 1, var )

d = ShuffleMatrix( x@dat)

t2 =apply( d, 1, var )

expect_true( !( all.equal(t1, t2) == TRUE),label="creating different results every time" )
