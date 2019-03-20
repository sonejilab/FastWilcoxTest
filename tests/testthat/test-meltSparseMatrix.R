context( 'meltSparseMatrix' )

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

meltedR = NULL;

system.time( {
	rownames(dat) = 1:nrow(dat)
	colnames(dat) = 1:ncol(dat)
	meltedR = reshape2::melt( dat )
	bad = which(meltedR[,3] == 0 )
	if ( length(bad) > 0)
		meltedR = meltedR[-bad,]
	
	}
)

system.time( { meltedCPP = meltSparseMatrix (x@dat) })

colnames(meltedR) = colnames( meltedCPP)

expect_equal( meltedR[,1], meltedCPP[,1] )
expect_equal( meltedR[,2], meltedCPP[,2] )
expect_equal( meltedR[,3], meltedCPP[,3] )

## 7.938 (R) vs. 0.429 (c++) using ncol = 10000 and nrow=9000