prefix='./'
load(file.path( prefix, 'data', 'pred.RData') )
load(file.path( prefix, 'data', 'inbag.RData') )

n = nrow(pred)
prox = matrix(0, n, n )
prox2 =  extract_proximity_oob ( pred, prox, inbag)

if ( ! file.exists( file.path( prefix, 'data', 'prox.RData') ) ) {
	for (i in 1:n) {
		for (j in 1:n) {
      		# Use only trees where both obs are OOB
      		tree_idx = which(inbag[i, ] == 0 & inbag[j, ] == 0)
      		if ( length(tree_idx) >0 ){
      			prox[i, j] = sum(pred[i, tree_idx] == pred[j, tree_idx]) / length(tree_idx)
   			}else {
   				prox[i, j] = 0
   			}
		}
	}
	save( prox, file = file.path( prefix, 'data', 'prox.RData'))
}

load(file.path( prefix, 'data', 'prox.RData') )

expect_equal( prox, prox2, label="proximity correct")
