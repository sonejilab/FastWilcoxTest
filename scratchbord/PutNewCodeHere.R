
Rstats <- function( X, interest, backgound,  logFCcut = 1.0 ){
	logRF = unlist(apply( X, 1, function( x ) { logFC( x[interest], x[backgound] ) } ) )
	OK = which( logFC > logFCcut )
	p.values = unlist(lapply( OK, function ( i ) {
			wilcox.test( X[interst, i], X[backgound,i], 'greater')$p.value
	} ))
	rbind( colID = OK, logFC = logFC[OK], p.value = p.values)
}
