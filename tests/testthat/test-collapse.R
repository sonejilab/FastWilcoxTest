
context( 'Cpp collapse sparse matrix')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

## nGenes for samples:
expect_equal( as.vector(apply( x@dat,2, function(d) { length(which(d > 0 )) } )), ColNotZero(x@dat) )

#### normal add - should be ~7x faster ( 11.829 / 1.648 )

ids = rep( 1:11, floor(ncol/11))
ids = c(ids, sample(1:11,ncol%%11))
ids = sort(ids)

expect_equal( as.vector(table(ids)), c(9,10,9,9,9,9,9,9,9,9,9), label="ids correct")

system.time({red = as.data.frame(collapse( x@dat, ids, 1 ))}) ## normal addition
colnames( red ) = 1:11;
rownames(red) = make.names(rownames(x@dat))

## now get the expected in R and make sure this is lower ;-)

sumUp <- function( x, ids ) {
	unlist(lapply( split( x, ids), function(d) { sum(d) } ))
}

system.time({ blue = t({ blue = data.frame(apply(  x@dat, 1, sumUp, ids ) )} ) } )

expect_equal(as.matrix(red), as.matrix(blue))

#### normal log add - should be ~6x faster ( 13.993 /2.335 )
system.time({green = as.data.frame(collapse( x@dat, ids, 2 ))}) ## normal addition
colnames( green ) = 1:11;
rownames(green) = make.names(rownames(x@dat))

sumUp <- function( x, ids ) {
	unlist(lapply( split( x, ids), function(d) { mean(d) } ))
}

system.time({ black = t({ black = data.frame(apply(  x@dat, 1, sumUp, ids ) )} ) } )

expect_equal(as.matrix(green), as.matrix(black))




