## R -d "valgrind --dsymutil=yes" 

devtools::load_all()
load('X.RData') # a big sparse matrix

Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) 

Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1)

## you will likely not reach this ;-)

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

message( 'Not broken after 4 runs - remarkable')

