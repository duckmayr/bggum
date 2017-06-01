load('newExData.RData')
library(ggumR)

chain1 <- ggumMCMC(responseMatrix, rep(4, 20), 50000)
cat('Chain 1 complete!\n')
chain2 <- ggumMCMC(responseMatrix, rep(4, 20), 50000)
cat('Chain 2 complete!\n')
chain3 <- ggumMCMC(responseMatrix, rep(4, 20), 50000)
cat('Chain 3 complete!\n')
save(chain1, chain2, chain3, file='testChains.RData')

tau0s <- seq(from=1041, to=1120, by=4)
tau1s <- tau0s + 1
tau2s <- tau1s + 1
tau3s <- tau2s + 1
inds <- c(seq(from=1, to=1000, by=10), seq(from=1001, to=50000, by=25))

pdf('thetaTrace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1:20, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Theta', x), xlab='Iteration', ylab='Value',
       ylim=c(-4, 4))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=thetas[x], lty=2)
}))
dev.off()

pdf('alphaTrace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1001:1020, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Alpha', x-1000), xlab='Iteration', ylab='Value',
       ylim=c(0, 4))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=alphas[x-1000], lty=2)
}))
dev.off()

pdf('deltaTrace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1021:1040, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Delta', x-1020), xlab='Iteration', ylab='Value',
       ylim=c(-5, 5))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=deltas[x-1020], lty=2)
}))
dev.off()

pdf('tau1Trace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau1s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Tau1', which(tau1s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-6, 6))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=taus[[which(tau1s == x)]][2], lty=2)
}))
dev.off()

pdf('tau2Trace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau2s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Tau2', which(tau2s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-6, 6))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=taus[[which(tau2s == x)]][3], lty=2)
}))
dev.off()

pdf('tau3Trace.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau3s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.5),
       main=paste('Tau3', which(tau3s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-6, 6))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.5))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.5))
  abline(h=taus[[which(tau3s == x)]][4], lty=2)
}))
dev.off()

