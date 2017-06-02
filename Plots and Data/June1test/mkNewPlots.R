load('June1chains.RData')
load('June1sims.RData')

tau0s <- seq(from=1041, to=1120, by=4)
tau1s <- tau0s + 1
tau2s <- tau1s + 1
tau3s <- tau2s + 1
inds <- c(seq(from=1, to=1000, by=10), seq(from=1001, to=50000, by=25))

pdf('thetaTraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1:20, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Theta', x), xlab='Iteration', ylab='Value',
       ylim=c(-3, 3))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=thetas[x], lty=2)
}))
dev.off()

pdf('thetaTrace2June1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(21:40, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Theta', x), xlab='Iteration', ylab='Value',
       ylim=c(-3, 3))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=thetas[x], lty=2)
}))
dev.off()

pdf('alphaTraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1001:1020, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Alpha', x-1000), xlab='Iteration', ylab='Value',
       ylim=c(0, 3))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=alphas[x-1000], lty=2)
}))
dev.off()

pdf('deltaTraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1021:1040, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Delta', x-1020), xlab='Iteration', ylab='Value',
       ylim=c(-2, 2))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=deltas[x-1020], lty=2)
}))
dev.off()

pdf('tau1TraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau1s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Tau1', which(tau1s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2, 2))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=taus[[which(tau1s == x)]][2], lty=2)
}))
dev.off()

pdf('tau2TraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau2s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Tau2', which(tau2s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2, 2))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=taus[[which(tau2s == x)]][3], lty=2)
}))
dev.off()

pdf('tau3TraceJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau3s, function(x){
  plot(inds, chain1[inds, x], type='l', col=rgb(0, 0, 1, alpha=0.4),
       main=paste('Tau3', which(tau3s == x)),
       xlab='Iteration', ylab='Value', ylim=c(-2, 2))
  lines(inds, chain2[inds, x], col=rgb(1, 0, 0, alpha=0.4))
  lines(inds, chain3[inds, x], col=rgb(0, 1, 0, alpha=0.4))
  abline(h=taus[[which(tau3s == x)]][4], lty=2)
}))
dev.off()

pdf('thetaDensityJune1.pdf', width=9, height=26)
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

pdf('thetaDensity2June1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(21:40, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Theta', x), xlab='Value', ylab='Density',
       xlim=c(-3, 3))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=thetas[x], lty=2)
}))
dev.off()

pdf('alphaDensityJune1.pdf', width=9, height=26)
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

pdf('deltaDensityJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(1021:1040, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Delta', x-1020), xlab='Value', ylab='Density',
       xlim=c(-2, 2))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=deltas[x-1020], lty=2)
}))
dev.off()

pdf('tau1DensityJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau1s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau1', which(tau1s == x)),
       xlab='Value', ylab='Density', xlim=c(-2, 2))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau1s == x)]][2], lty=2)
}))
dev.off()

pdf('tau2DensityJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau2s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau2', which(tau2s == x)),
       xlab='Value', ylab='Density', xlim=c(-2, 2))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau2s == x)]][3], lty=2)
}))
dev.off()

pdf('tau3DensityJune1.pdf', width=9, height=26)
layout(matrix(1:20, ncol=2))
invisible(sapply(tau3s, function(x){
  plot(density(chain1[ , x]), col='blue',
       main=paste('Tau3', which(tau3s == x)),
       xlab='Value', ylab='Density', xlim=c(-2, 2))
  lines(density(chain2[ , x]), col='red')
  lines(density(chain3[ , x]), col='green')
  abline(v=taus[[which(tau3s == x)]][4], lty=2)
}))
dev.off()

