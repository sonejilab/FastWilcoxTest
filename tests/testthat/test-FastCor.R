

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

context( 'roll sum')
m = Matrix::Matrix( 1:10, nrow=1, sparse=T)
all.equal( rollSum(m,3 ), matrix(c(6,9,12,15,18,21,24,27), nrow=1))

rollSumR <- function( x, n) {
	unlist(lapply( (n+1):(length(x)+1) , function(a) { sum(x[(a-n):(a-1)] )} ))
}

all.equal( rollSum(m,3 ), matrix(rollSumR(m[1,],3), nrow=1))



n = 10
system.time({ rolled = rollSum( x@dat, n ) })

#colnames(rolled) = colnames(x@dat)[(n+1):ncol(x@dat)]
#rownames(rolled) = rownames(x@dat)

system.time({ cmp = t(apply( x@dat, 1, rollSumR, n) ) })
#colnames(cmp) = colnames(x@dat)[(n+1):ncol(x@dat)]
rownames(cmp) =NULL
expect_equal( as.vector(rolled), as.vector(cmp))

system.time( { stats = CorNormalMatrix(  rolled, 1:nrow(rolled) ) })

system.time({ exp =  cor( rolled, 1:nrow(rolled)) })

expect_equal( stats, as.vector(exp))

start=1:100
s = sample(20:40, 12)
x@dat[1:20,1:20] = x@dat[1:20,1:20] * 2

res = t(rollSumStart( Matrix::t(x@dat), 10, start))
