
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
all.equal( as.vector(apply( x@dat,2, function(d) { length(which(d > 0 )) } )), toColNotZero(x@dat) )



#### normal add - should be ~7x faster ( 11.829 / 1.648 )

ids = rep( 1:10, ceiling(ncol/10))

system.time({red = as.data.frame(collapse( x@dat, ids, 1 ))}) ## normal addition
colnames( red ) = 1:10;
rownames(red) = make.names(rownames(x@dat))

## now get the expected in R and make sure this is lower ;-)

sumUp <- function( x, ids ) {
	unlist(lapply( split( x, ids), function(d) { sum(d) } ))
}

system.time({ blue = t({ blue = data.frame(apply(  x@dat, 1, sumUp, ids ) )} ) } )

all.equal(as.matrix(red), as.matrix(blue))


#### normal log add - should be ~6x faster ( 13.993 /2.335 )

skip('broken')
system.time({red = as.data.frame(collapse( x@dat, ids, 0 ))}) ## log addition
colnames( red ) = 1:10;
rownames(red) = make.names(rownames(x@dat))

sumUp <- function( x, ids ) {
	log(unlist(lapply( split( x, ids), function(d) { sum(exp(d[which(d > 0 )])) } )))
}

system.time({ blue = t({ blue = data.frame(apply(  x@dat, 1, sumUp, ids ) )} ) } )

all.equal(as.matrix(red), as.matrix(blue), 1e-7)


#### normal std::expm1 add - should be ~6x faster ( 14.011 / 2.216 )

system.time({red = as.data.frame(collapse( x@dat, ids, 2 ))}) ## log addition
colnames( red ) = 1:10;
rownames(red) = make.names(rownames(x@dat))

sumUp <- function( x, ids ) {
	log(unlist(lapply( split( x, ids), function(d) { sum(exp(d[which(d > 0 )]) -1) } )))
}

system.time({ blue = t({ blue = data.frame(apply(  x@dat, 1, sumUp, ids ) )} ) } )

all.equal(as.matrix(red), as.matrix(blue), 1e-7)


