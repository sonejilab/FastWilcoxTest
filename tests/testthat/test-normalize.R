
context( 'Cpp collapse sparse matrix')

set.seed(1)
ncol = 1000
nrow=9000
dat = matrix(round(rnorm(ncol*nrow,mean = 5, sd = 7)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

A <- as_FastWilcoxTest( dat )

rm(dat)

to = min(Matrix::colSums(A@dat)) 

ret <- NormalizeCells(A@dat, to, display_progress = FALSE) ; 

#browser()
expect_equal( length(which(ret@i == -1)) , 0)

value =apply( ret, 2, function(x) sum(x[which(x > 0 )]) )

expect_equal( value, rep(to, ncol) )


## now I need to check, that the values actually make sens!
d= lapply ( 1:100, function( i ) {
	raw = as.vector(A@dat[,i])
	norm = ret[,i]
	lost = which(norm == -1 )
	norm_double = raw / sum(raw) * to
	## have all values that have been lost been below 0 after the norm?
	expect_equal( norm_double[lost] , rep( 0, length(lost)),1)
	## are all other values in the range +-1 of the norm data?
	expect_equal(norm[-lost], norm_double[-lost], 1 )
	NULL;
		} 
)

## use inbult data

ret <- NormalizeCells(x@dat, 500, display_progress = FALSE) ; 

sums = apply( ret ,2, function(x) sum(x[which(x > 0 )] ) )
exp = rep( 500, 100 )
exp[c( 8, 12, 13, 19, 22, 29, 41, 42, 44, 46, 51, 55, 56, 59, 62, 66, 70, 73, 79, 81, 85, 94, 95, 97 )] = 0
expect_equal( sums , exp,  0)



