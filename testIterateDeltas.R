source('~/ggum/ggumR/R/getPrior.R')
source('~/ggum/ggumR/R/probability.R')
load('~/ggum/exampleData.RData')
source('~/ggum/ggumR/R/Beta_ab.R')
source('~/ggum/ggumR/R/llCol.R')
source('~/ggum/ggumR/R/iterateDeltas.R')
source('~/ggum/ggumR/R/proposer.R')

# Function to give the rate of proposal acceptance for one parameter:
acceptanceRate <- function(chainVec){
  result <- 0
  for ( i in 2:length(chainVec) ) {
    if ( chainVec[i] != chainVec[i-1] ) {
      result <- result + 1
    }
  }
  return(result/length(chainVec))
}

# Do 5000 iterations of the algorithm for deltas given true values for all the
# other parameters:
testMat <- iterateDeltas(thetas, alphas, taus, 5000, responseMatrix[1:200, ], 4)
# Find the acceptance rate:
apply(testMat, 2, acceptanceRate)
# How are the estimates?
chainLength <- nrow(testMat)
ests <- apply(testMat[(chainLength/2):chainLength, ], 2, acceptanceRate)
mean(abs(ests - deltas)) # Not great
# Save the results for future examination:
save(testMat, file='testMatIterateDeltas.RData')

# Let's look at densities and trace plots:
pdf('deltaDensities.pdf', width=9, height=13)
layout(matrix(1:10, ncol=2))
invisible(sapply(1:10, function(x){
  plot(density(testMat[(chainLength/2):chainLength, x]),
       main=paste('Delta', x))
  abline(v=deltas[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Density', 'True Value'))
}))
dev.off()

pdf('deltaTracePlots.pdf', width=9, height=13)
layout(matrix(1:10, ncol=2))
invisible(sapply(1:10, function(x){
  Title <- paste('Trace to Iteration 1000', paste0('For Delta ', x), sep='\n')
  plot(testMat[1:1000, x], type='l', main=Title, xlab='Iteration', ylab='Value')
  abline(h=deltas[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Trace', 'True Value'))
}))
dev.off()