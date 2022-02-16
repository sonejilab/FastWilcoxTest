set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

#warning("Still OK 1")
d1 = ShuffleMatrix( x@dat)

#warning("Still OK 2")
t1 =apply( d1, 1, var )
#Sys.sleep(1)
#warning("Still OK 3")
d2 = ShuffleMatrix( x@dat)

#warning("Still OK 4")
t2 =apply( d2, 1, var )

expect_true( ! all ( t2 == t1) ,label="creating different results every time" )
