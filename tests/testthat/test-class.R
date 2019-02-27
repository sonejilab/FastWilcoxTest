context( 'Class usage')

set.seed(1)
dat = matrix(round(rnorm(9000,mean = 1, sd = 1)),ncol=100)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:100)
rownames(dat) <- paste( 'gene', 1:90)

x <- as_RcppTestArea( dat )

expect_equal( class( x )[1], c('RcppTestArea') )

A = 1:6
B = 2:8
E = log( mean(A) / mean(B) )

a= log(A)
b = log(B)

expect_equal( logFC ( a, b ), E)

system.time({ Cpp = StatTest(Matrix::t(x@dat), 1:20, 30:80, .1) })
system.time({ R = Rstats(Matrix::t(x@dat), 1:20, 30:80, .1) })
