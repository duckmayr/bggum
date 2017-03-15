source('probability.R')
source('probMatrix.R')
source('stableProd.R')
source('likelihood.R')
source('Beta_ab.R')
source('getPrior.R')
source('ratio.R')
source('proposer.R')
source('decider.R')
source('iteration.R')

## To make some fake data using the most likely response given parameters:
likelyResponse <- function(theta, alpha, delta, tau){
  K <- length(tau)
  return(which.max(sapply(1:K, probability, K, theta, alpha, delta, tau)))
}

likelyResponseVector <- function(theta, alpha, delta, tau){
  return(sapply(1:length(tau), function(x){
    likelyResponse(theta, alpha[x], delta[x], tau[[x]])}))
}

alphaDistParams <- list(shape1=1.5, shape2=1.5, a=.25, b=4)
deltaDistParams <- list(shape1=2, shape2=2, a=-5, b=5)
tauDistParams <-list(shape1=2, shape2=2, a=-6, b=6)
thetas <- rnorm(10)
alphas <- rBeta_ab(3, params=alphaDistParams)
deltas <- rBeta_ab(3, params=deltaDistParams)
taus <- list(c(0, rBeta_ab(2, params=tauDistParams)),
             c(0, rBeta_ab(3, params=tauDistParams)),
             c(0, rBeta_ab(4, params=tauDistParams)))
responseMat <- matrix(sapply(thetas, likelyResponseVector, alphas, deltas,
                             taus), nrow=10, byrow=TRUE)

## Now demonstrate the functions:
resultList <- iteration(thetas, alphas, deltas, taus, responseMat)
for (i in 1:10) {
  print(c(thetas[i], resultList$thetas[i]))
}
for (i in 1:3) {
  print(c(alphas[i], resultList$alphas[i]))
}
for (i in 1:3) {
  print(c(deltas[i], resultList$deltas[i]))
}
for (i in 1:3) {
  for (j in 1:length(taus[[i]])) {
    print(c(taus[[i]][j], resultList$taus[[i]][j]))
  }
}
