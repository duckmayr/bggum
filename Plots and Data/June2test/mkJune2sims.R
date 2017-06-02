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
alphas <- rtruncnorm(20, 1, 0.5, 0.25, 4)
deltas <- c(rtruncnorm(14, 0, 0.5, -5, 5), rep(3, 3), rep(-3, 3))
taus <- as.list(sapply(1:20, function(x){
  c(0, rtruncnorm(3, 0, 0.5, -6, 6))
}, simplify=FALSE))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)
apply(responseMatrix, 2, table)
# Item 20 only had responses at options 1 and 2
taus[[20]] <- c(0, rtruncnorm(3, 0, 0.5, -6, 6))
responseMatrix[ , 20] <- sapply(1:1000, function(x){
  guessResponse(thetas[x], alphas[20], deltas[20], taus[[20]])
})
apply(responseMatrix, 2, table)
# That fixed it, though the items with extreme delta values
# (that should therefore be montonic)
# had sparser responses for some options
save(thetas, alphas, deltas, taus, responseMatrix,
     file='June2sims.RData')
