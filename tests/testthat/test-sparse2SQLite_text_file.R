context( 'sparse2SQLite_text_file' )

prefix = './'
#prefix = './tests/testthat/'

set.seed(1)
ncol = 1000
nrow=9000
dat = matrix(round(rnorm(ncol*nrow,mean = 5, sd = 7)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

A <- as_FastWilcoxTest( dat )

ofile = file.path( prefix, 'data', 'Output', 'cpp_melted_table.txt')
if ( file.exists(ofile)){
	unlink(ofile)
}
start_time <- Sys.time()
sparse2SQLite_text_file( A@dat, ofile )
end_time <- Sys.time()

print ( paste("load cellexalObj:",difftime(end_time, start_time,  units = "secs")[[1]], "sec"))

expect_true( file.exists( ofile), label="c++ has created the oputfile" )

dat = read.delim(header=F, file= ofile, sep= " " )

melted = meltSparseMatrix( A@dat )

colnames(dat) = colnames(melted)

expect_equal( dat, melted, label="c++ file is the same as melted in R" )