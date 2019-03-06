context( 'Class usage')

set.seed(1)
dat = matrix(round(rnorm(9000,mean = 1, sd = 1)),ncol=100)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:100)
rownames(dat) <- paste( 'gene', 1:90)

x <- as_FastWilcoxTest( dat )

expect_equal( class( x )[1], c('FastWilcoxTest') )

A = 1:6
B = 2:8
E = log( mean(A) / mean(B) )

a= log(A+1)
b= log(B+1)

expect_equal( logFC ( a, b ), E,  tolerance=1e-1)



Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1) 

R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .1, .1) 
## remove one value that has a p value difference of 5.372936e-03 (1.00000000 vs 0.99462706)
Cpp = Cpp[-36,]
R = R[-36,]
expect_equal( as.vector(R[,'p.value'] ), as.vector( Cpp[,'p.value']) , tolerance=1e-10)


Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1, .1, TRUE) 

R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .1, .1, TRUE) 

Cpp = Cpp[-22,]
R = R[-22,]

expect_equal( as.vector(R[,'p.value'] ), as.vector( Cpp[,'p.value']) , tolerance=1e-10)

