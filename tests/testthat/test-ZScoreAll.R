context( 'Cpp z.score')

set.seed(1)
ncol = 100
nrow=900
dat = matrix(round(rnorm(ncol*nrow,mean = 3, sd = 5)),ncol=ncol)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)

x <- as_FastWilcoxTest( dat )

rm(dat)

bad = sample( 1:length( x@dat@x), round( length( x@dat@x) / 100) )  
x@dat@x[ bad ] = -1 # simulate my normalization result

system.time({zscored = ZScoreAll(x@dat)})

system.time({zscoredR = RzscoreAll(x)})

rownames(zscored) = rownames(x@dat)
colnames(zscored) = colnames(x@dat)

expect_equal( zscored, zscoredR, label="faster zscore method works correctly" )