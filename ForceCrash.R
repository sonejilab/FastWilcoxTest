devtools::load_all()
load('X.RData') # a big sparse matrix

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

system.time({ Cpp = StatTest(Matrix::t(X), 1:20, 30:80, .1) })

message( 'Not broken after 4 runs - remarkable')

