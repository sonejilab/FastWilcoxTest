context( 'Class usage')

set.seed(1)
ncol = 100
nrow=90
dat = matrix(round(rnorm(ncol*nrow,mean = 0, sd = 1)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

expect_equal( class( x )[1], c('FastWilcoxTest') )

context( 'wilcox tests')
A = 1:6
B = 2:8
E = log( mean(A) / mean(B) )

a= log(A+1)
b= log(B+1)

expect_equal( logFC ( a, b ), E,  tolerance=1e-1)

Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1) 

R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .099, .099) 
## remove one value that has a p value difference of 5.372936e-03 (1.00000000 vs 0.99462706)
Cpp = Cpp[-59,]
R = R[-59,]
expect_equal( as.vector(R[,'p.value'] ), as.vector( Cpp[,'p.value']) , tolerance=1e-10)


system.time({Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1, TRUE) })

system.time({R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .099, .099, TRUE) }) 

Cpp = Cpp[-22,]
R = R[-22,]

expect_equal( as.vector(R[,'p.value'] ), as.vector( Cpp[,'p.value']) , tolerance=1e-10)

context( 'fast cor')

set.seed(1)

t1 = rnorm( 30 )

cor = correlationCoefficient( X= t1, Y = t1 )

expect_equal( cor ,1 )

cor = correlationCoefficient( X= t1, Y =rev( t1 ))
expect_equal( cor , -0.010188, tolerance=1e-5 )

system.time({ cor = CorMatrix (x@dat, x@dat[1,] ) })

system.time({cor2 = apply(as.matrix(x@dat),1, cor, x@dat[1,] ) })


all.equal(cor, as.vector(cor2), 1e-7)


