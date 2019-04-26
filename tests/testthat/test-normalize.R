
context( 'Cpp collapse sparse matrix')

set.seed(1)
ncol = 1000
nrow=9000
dat = matrix(round(rnorm(ncol*nrow,mean = 5, sd = 7)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

to = min(Matrix::colSums(x@dat)) 

ret <- NormalizeCells(x@dat, to, display_progress = FALSE) ; 

expect_equal( length(which(ret@i == -1)) , 0)

value =apply( ret, 2, function(x) sum(x[which(x > 0 )]) )

expect_equal( value, rep(to, ncol) )


## now I need to check, that the values actually make sens!
d= lapply ( 1:100, function( i ) {
	raw = as.vector(x@dat[,i])
	norm = ret[,i]
	lost = which(norm == -1 )
	norm_double = raw / sum(raw) * to
	## have all values that have been lost been below 0 after the norm?
	expect_equal( norm_double[lost] , rep( 0, length(lost)),1)
	## are all other values in the range +-1 of the norm data?
	expect_equal(norm[-lost], norm_double[-lost], 1 )
	NULL;
		} )


context( 'Sample scale Lg norm')


dat = matrix(round(rnorm(ncol*nrow,mean = 5, sd = 7)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

ret <- NormalizeSamples(x@dat, rep( 100, ncol) , display_progress = FALSE) ;
colnames(ret) = colnames(x@dat)
rownames(ret) = rownames(x@dat)

expect_true( all.equal( as.matrix(x@dat) / 100, as.matrix(ret)) )



