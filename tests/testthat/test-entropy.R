context( 'entropy')

skip('useless and depricated')

# entR <- function (d) {
#     if ( length(d) < 5) {
#       return( 0 )
#     }

#   a = table(d)
#   ret = 0
#    if ( length(a) > 1 ) {
#       p = a/sum(a)
#       ret = -sum(p*log(p))
#    }
#  ret
# }

# set.seed ( 1 )
# dat = sample( 1:10, 1e+8, replace=T );

# start_time <- Sys.time()
# entro = entropy( dat )
# end_time <- Sys.time()
# print ( paste("c++ version:",difftime(end_time, start_time,  units = "secs")[[1]], "sec"))


# start_time <- Sys.time()
# entro = entroR = entR( dat )
# end_time <- Sys.time()
# print ( paste("R version",difftime(end_time, start_time,  units = "secs")[[1]], "sec"))


# expect_equal( entro, entroR)

# x = matrix(runif( 3000 ), ncol=3)

# X = c( min(x[,1]), max(x[,1]) )
# Y = c( min(x[,2]), max(x[,2]) )
# Z = c( min(x[,3]), max(x[,3]) )
# dist = FastWilcoxTest::euclidian_distances3d( X, Y ,Z )

# expect_equal( dist[2], 1.729509, 1e-4)

# radii = c( dist[2]/100, dist[2]/10, dist[2] )

# system.time({entropyT = SphericEntropy( x[,1], x[,2], x[,3], dat[1:1000], radii)})

# SphericEntropyR <- function (x1,x2,x3,gvect, n ) {
# 	closest = NULL
#     for( res in  lapply( 1:length(gvect), function( id ) {
#         d = FastWilcoxTest::eDist3d( x1, x2, x3, id -1 )
#         order(d)[1:n]
#         m = matrix(
#             unlist( lapply( n , function( N ) {
#                 OK = which(d < N)
#                 if ( length(OK) == 0 ){q
#                     return(0)
#                 }
#             c( FastWilcoxTest::entropy(  gvect[ OK ]), length(OK) )
#             } ))
#         , nrow=2)
#         rownames(m) = paste( 'cell', id, c('entropy', 'selectedCells'))
#         m
#      })) {
#         closest = rbind( closest, res )
#     }
#     closest
# }

# entropyR = SphericEntropyR( x[,1], x[,2], x[,3], dat[1:1000], radii)

# colnames(entropyR) = colnames(entropyT)

# all.equal( entropyR, entropyT )
