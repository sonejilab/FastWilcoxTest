
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

context( 'roll sum Area')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rolled = rollAreaSum ( x@dat, location=1:90, size=10, funcID=1 )
#browser()
N = c(rep( 9, 80) )
rolledR = apply(x@dat, 2, function(Z){
	ret = vector('numeric', length(N))
	for ( i in 1:length(N)) { 
		ret[i] = sum(Z[i:(i+N[i])])
	}
	ret
	} 
) 
expect_true(all(rolled== rolledR ))