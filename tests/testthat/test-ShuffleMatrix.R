set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )
xv = apply( x@dat[,1:50], 2, var )
xn = FastWilcoxTest::ColNotZero( x@dat[,1:50])
rm(dat)

#warning("Still OK 1")
d1 = ShuffleMatrix( x@dat)

expect_equal( xn, FastWilcoxTest::ColNotZero( d1) ,label="same density first 50 cols" )

tr1 =apply( d1, 2, var )
expect_equal( as.vector(xv), tr1 ,label="same variance (cells) first 50 cols" )

t1 =apply( d1, 1, var )
expect_true( ! all (apply( x@dat[,1:50], 1, var ) ==  t1) ,label="different variance (cells) first 50 cols" )

#Sys.sleep(1)
#warning("Still OK 3")
d2 = ShuffleMatrix( x@dat)

#warning("Still OK 4")
t2 =apply( d2, 1, var )

expect_true( ! all ( t2 == t1) ,label="creating different results every time" )

prefix = "."
path = file.path( prefix, 'data', 'problematic.RData' )


load( path )
dat = Matrix::drop0( x@dat)
d3 = ShuffleMatrix( dat)

tr3 =apply( d3, 2, var )
expect_equal( as.vector(apply(dat,2, var)), tr3 ,label="same variance (cells) first 50 cols" )

