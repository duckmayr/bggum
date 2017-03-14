source('probability.R') # source in the probability function

# The following are the parameters for the example in de la Torre (2006):
a <- 2; d <- 0; t <- c(0, -1, -0.7, -0.4); th <- seq(from=-3, to=3, by=0.1)
# So we store in an object the probability an individual would choose each
# of the item's four options at various thetas:
probs <- t(sapply(th, function(x){
  sapply(1:length(t), probability, length(t), x, a, d, t)
}))
# And plot the option response functions:
plot(NULL, ylim=c(0, 1), xlim=c(-3, 3), ylab='Prob. of Positive Response',
     xlab=expression(theta), main='GGUM Option Response Functions')
plotColors <- c('grey40', 'grey50', 'grey20', 'grey30')
plotSymbols <- c(18, 20, 17, 15)
indices <- c(TRUE, FALSE)
matlines(th, probs, col=plotColors, lty=1)
matpoints(th[indices], probs[indices, ], col=plotColors, pch=plotSymbols)
legend('right', legend=c('Opt. 1', 'Opt. 2', 'Opt. 3', 'Opt. 4'),
       pch=plotSymbols, col=plotColors)