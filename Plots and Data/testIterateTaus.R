source('~/ggum/ggumR/R/getPrior.R')
source('~/ggum/ggumR/R/probability.R')
load('~/ggum/exampleData.RData')
source('~/ggum/ggumR/R/Beta_ab.R')
source('~/ggum/ggumR/R/llCol.R')
source('~/ggum/ggumR/R/iterateTaus.R')
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

# Do 5000 iterations of the algorithm for taus given true values for all the
# other parameters:
testMat <- iterateTaus(thetas, alphas, deltas, 5000, responseMatrix[1:200, ], 4)
# Find the acceptance rate:
apply(testMat, 2, acceptanceRate)
# How are the estimates?
chainLength <- nrow(testMat)
ests <- apply(testMat[(chainLength/2):chainLength, ], 2, mean)
mean(abs(ests - unlist(taus))) # Not too bad
# Save the results for future examination:
save(testMat, file='testMatIterateTaus.RData')

# Let's look at densities and trace plots:
pdf('tauDensities.pdf', width=18, height=25)
layout(matrix(1:40, ncol=4, byrow=TRUE))
invisible(sapply(1:40, function(x){
  plot(density(testMat[(chainLength/2):chainLength, x]),
       main=paste0('Tau ',((x - 1) %/% 4) + 1, '-', ((x - 1) %% 4) + 1))
  abline(v=unlist(taus)[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Density', 'True Value'))
}))
dev.off()

pdf('tauTracePlots.pdf', width=18, height=25)
layout(matrix(1:40, ncol=4, byrow=TRUE))
invisible(sapply(1:40, function(x){
  Title <- paste('Trace to Iteration 1500',
                 paste0('For Tau ', ((x - 1) %/% 4) + 1, '-',
                        ((x - 1) %% 4) + 1), sep='\n')
  plot(testMat[1:1500, x], type='l', main=Title, xlab='Iteration', ylab='Value')
  abline(h=unlist(taus)[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Trace', 'True Value'))
}))
dev.off()