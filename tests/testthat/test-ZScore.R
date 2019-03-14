context( 'Cpp z.score')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 3, sd = 5)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

zscored = ZScore(x@dat)

expect_equal( dim(x@dat), dim(zscored) )

Rinfo = cbind(  apply(x@dat,1,function(d) { mean(d[which(d>0)] ) } ),  
		apply(x@dat,1,function(d) { x= d[which(d>0)]; sum( (x - mean(x)) * (x-mean(x)) ) } ),  
		apply(x@dat,1,function(d) { sd(d[which(d>0)] ) } )
)

CPPinfo = MEAN_STD( x@dat )

expect_equal( as.vector(Rinfo[,1]), as.vector(CPPinfo[,1]) )
expect_equal( as.vector(Rinfo[,2]), as.vector(CPPinfo[,2]) )
expect_equal( as.vector(Rinfo[,3]), as.vector(CPPinfo[,3]) )


## R z.score...
z.score = function ( x ){
	
	i = i+1
	n <- which(x <= 0)
	dropped = which(x == -1)
	if ( length(x) - length(n) > 1 ){
		if (length(n) == 0 ){
			x <-  10 + scale(as.vector(t(x)))
		}
		else {
			x[-n] <- 10 + scale(as.vector(t(x[-n])))
		}
		
	}
	else {
	}
	x
}

i=0
zscoredR =  t(apply( x@dat,1, z.score))

expect_equal( dim(zscoredR), dim(zscored) )

expect_equal( Matrix::Matrix(zscoredR, sparse=T)@x , zscored@x, 6e-4 )
