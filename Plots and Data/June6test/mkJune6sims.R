library(ggumR)

guessResponse <- function(theta, alpha, delta, tau){
  K <- length(tau)
  probs <- ggumProbability(1:K, theta, alpha, delta, tau)
  randomNumber <- runif(1)
  for (i in (K-1):1) {
    if (randomNumber > sum(probs[1:i])) {
      return(i + 1)
    }
  }
  return(1)
}

mkResponseMatrix <- function(thetas, alphas, deltas, taus){
  n <- length(thetas)
  m <- length(alphas)
  responseMatrix <- matrix(0, nrow=n, ncol=m)
  for (i in 1:n) {
    for (j in 1:m) {
      responseMatrix[i, j] <- guessResponse(thetas[i], alphas[j], deltas[j],
                                            taus[[j]])
    }
  }
  return(responseMatrix)
}

thetas <- rnorm(1000)
alphas <- rtruncnorm(20, 1, 1, 0.25, 4)
deltas <- c(rtruncnorm(14, 0, 1, -5, 5), seq(from=2.75, to=3.25, by=0.25),
            seq(from=-2.75, to=-3.25, by=-0.25))
taus <- as.list(sapply(1:20, function(x){
  c(0, runif(3, -2, 0))
}, simplify=FALSE))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)
apply(responseMatrix, 2, table)
save(thetas, alphas, deltas, taus, responseMatrix,
     file='June6sims.RData')
