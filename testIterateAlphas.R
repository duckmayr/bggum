source('~/ggum/ggumR/R/getPrior')
source('~/ggum/ggumR/R/probability.R')
load('~/ggum/exampleData.RData')
source('~/ggum/ggumR/R/Beta_ab.R')
source('~/ggum/ggumR/R/llCol.R')
source('~/ggum/ggumR/R/iterateAlphas.R')
source('~/ggum/ggumR/R/proposer.R')
acceptanceRate <- function(chainVec){
  result <- 0
  for ( i in 2:length(chainVec) ) {
    if ( chainVec[i] != chainVec[i-1] ) {
      result <- result + 1
    }
  }
  return(result/length(chainVec))
}
testMat <- iterateAlphas(thetas, deltas, taus, 5000, responseMatrix[1:200, ], 4)
apply(testMat, 2, acceptanceRate)
chainLength <- nrow(testMat)
ests <- apply(testMat[(chainLength/2):chainLength, ], 2, acceptanceRate)
mean(abs(ests - alphas))
save(testMat, file='testMatIterateAlphas.RData')

pdf('alphaDensities.pdf', width=9, height=13)
layout(matrix(1:10, ncol=2))
invisible(sapply(1:10, function(x){
  plot(density(testMat[(chainLength/2):chainLength, x]),
       main=paste('Alpha', x))
  abline(v=alphas[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Density', 'True Value'))
}))
dev.off()

pdf('alphaTracePlots.pdf', width=9, height=13)
layout(matrix(1:10, ncol=2))
invisible(sapply(1:10, function(x){
  Title <- paste('Trace to Iteration 1000', paste0('For Alpha ', x), sep='\n')
  plot(testMat[1:1000, x], type='l', main=Title, xlab='Iteration', ylab='Value')
  abline(h=alphas[x], lty=2)
  legend('topright', bty='n', lty=c(1, 2), legend=c('Trace', 'True Value'))
}))
dev.off()