context( 'LinLang tests')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)
group = sample( 1:5, ncol, replace = TRUE)

x <- as_FastWilcoxTest( dat )

rm(dat)


system.time( {res = LinLang(Matrix::t(x@dat), Grouping= group, nGroup=5, display_progress=TRUE )})


RLinLang <- function( x, Grouping, nGroup) {
	d = split( x, Grouping)
	d1 = data.frame(lapply( d, function(x) {
		c( length(x), length(which(x!= 0)), mean( x[which(x!= 0)]) )
	}))
	d1[3,which(is.na(d1[3,]))] = 0
	cmp = 1:nGroup
	A = cor(cmp,  as.vector(t(d1[2,]) / t(d1[1,])) )
	B = cor( cmp , as.vector(t(d1[3,])))
	C = cor( as.vector(t(d1[2,]) / t(d1[1,])), as.vector(t(d1[3,])))
	c( A, B, C, var( as.vector(t(d1[2,]) / t(d1[1,]))) , var( as.vector(t(d1[3,])))) 
}

system.time( {resR = t(data.frame(apply( x@dat, 1, RLinLang, Grouping= group, nGroup=5 )))})
rownames(resR) = NULL
colnames(resR) = colnames(res)
expect_equal( res , as.matrix(resR) ,1e-7)
