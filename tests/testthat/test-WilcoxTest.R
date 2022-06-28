context( 'wilcox tests')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

browser()
x <- as_FastWilcoxTest( dat )

rm(dat)

A = 1:6
B = 2:8
E = mean(A) / mean(B) 

expect_equal( FC ( B, A ), E)

Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1) 

fracExprIN = ColNotZero(  Matrix::t(x@dat[,1:20])) / 20
fracExprOUT = ColNotZero(  Matrix::t(x@dat[,30:80])) / length(30:80)

expect_true( length(which(! fracExprIN[Cpp[,1]] ==  Cpp[,3])) ==0, label="fracExprIN correct" )
expect_true( length(which(! fracExprOUT[Cpp[,1]] ==  Cpp[,4])) ==0, label="fracExprOUT correct" )

fracExprIN = apply ( x@dat[,1:20], 2, function(x) { length(which(x>0)) / 20 } )
R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .09, .09) 

expect_equal( as.vector(R[,'p.value'][-70] ), as.vector( Cpp[,'p.value'])[-70])

system.time({Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1, TRUE) })
system.time({R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .099, .099, TRUE) })

Cpp = Cpp[-70,]
R = R[-70,]

expect_equal( as.vector(R[,'p.value'] ), as.vector( Cpp[,'p.value']) )

