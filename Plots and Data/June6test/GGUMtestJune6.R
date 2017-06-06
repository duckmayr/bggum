load('June6sims.RData')
library(ggumR)

chain1 <- ggumMCMC(responseMatrix, rep(4, 20), 30000)
cat('Chain 1 complete!\n')
chain2 <- ggumMCMC(responseMatrix, rep(4, 20), 30000)
cat('Chain 2 complete!\n')
chain3 <- ggumMCMC(responseMatrix, rep(4, 20), 30000)
cat('Chain 3 complete!\n')
save(chain1, chain2, chain3, file='June6chains.RData')

tau0s <- seq(from=1041, to=1120, by=4)
tau1s <- tau0s + 1
tau2s <- tau1s + 1
tau3s <- tau2s + 1
inds <- c(seq(from=1, to=1000, by=10), seq(from=1001, to=50000, by=25))

pdf('thetaTraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1:20, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Theta', x), xlab='Iteration', ylab='Value',
       ylim=c(-3.75, 3.75))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=thetas[x], lty=2)
}))
dev.off()

pdf('alphaTraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1001:1020, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Alpha', x-1000), xlab='Iteration', ylab='Value',
       ylim=c(0, 4))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=alphas[x-1000], lty=2)
}))
dev.off()

pdf('deltaTraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1021:1040, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Delta', x-1020), xlab='Iteration', ylab='Value',
       ylim=c(-3.5, 3.5))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=deltas[x-1020], lty=2)
}))
dev.off()

pdf('tau1TraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau1s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Tau1', which(tau1s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2.25, 0.25))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=taus[[which(tau1s == x)]][2], lty=2)
}))
dev.off()

pdf('tau2TraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau2s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Tau2', which(tau2s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2.25, 0.25))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=taus[[which(tau2s == x)]][3], lty=2)
}))
dev.off()

pdf('tau3TraceJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau3s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.45),
       main=paste('Tau3', which(tau3s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2.25, 0.25))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.45))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.45))
  abline(h=taus[[which(tau3s == x)]][4], lty=2)
}))
dev.off()

pdf('thetaDensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1:20, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Theta', x), xlab='Value', ylab='Density',
       xlim=c(-3, 3))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=thetas[x], lty=2)
}))
dev.off()

pdf('alphaDensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1001:1020, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Alpha', x-1000), xlab='Value', ylab='Density',
       xlim=c(0, 3))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=alphas[x-1000], lty=2)
}))
dev.off()

pdf('deltaDensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1021:1040, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Delta', x-1020), xlab='Value', ylab='Density',
       xlim=c(-3.5, 3.5))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=deltas[x-1020], lty=2)
}))
dev.off()

pdf('tau1DensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau1s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau1', which(tau1s == x)),
       xlab='Value', ylab='Density', xlim=c(-2.25, 0.25))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau1s == x)]][2], lty=2)
}))
dev.off()

pdf('tau2DensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau2s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau2', which(tau2s == x)),
       xlab='Value', ylab='Density', xlim=c(-2.25, 0.25))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau2s == x)]][3], lty=2)
}))
dev.off()

pdf('tau3DensityJune6.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau3s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau3', which(tau3s == x)),
       xlab='Value', ylab='Density', xlim=c(-2.25, 0.25))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau3s == x)]][4], lty=2)
}))
dev.off()
