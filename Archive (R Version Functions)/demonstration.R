setwd('~/ggum')
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
guessResponse <- function(theta, alpha, delta, tau){
  K <- length(tau)
  probs <- probability(1:K, K, theta, alpha, delta, tau)
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

alphaDistParams <- list(shape1=1.5, shape2=1.5, a=.25, b=4)
deltaDistParams <- list(shape1=2, shape2=2, a=-5, b=5)
tauDistParams <-list(shape1=2, shape2=2, a=-6, b=6)
thetas <- rnorm(10)
alphas <- rBeta_ab(3, params=alphaDistParams)
deltas <- rBeta_ab(3, params=deltaDistParams)
taus <- list(c(0, rBeta_ab(2, params=tauDistParams)),
             c(0, rBeta_ab(3, params=tauDistParams)),
             c(0, rBeta_ab(4, params=tauDistParams)))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)

## Now demonstrate the functions:
ggum <- function(responseMatrix, Kvector, iterations){
  n <- nrow(responseMatrix)
  m <- ncol(responseMatrix)
  initialThetas <- rnorm(n)
  initialAlphas <- rBeta_ab(m, params=list(shape1=1.5, shape2=1.5, a=.25, b=4))
  initialDeltas <- rBeta_ab(m, params=list(shape1=2, shape2=2, a=-5, b=5))
  initialTaus <- as.list(sapply(1:m, function(x){
    c(0, rBeta_ab(Kvector[x]-1, params=list(shape1=2, shape2=2, a=-6, b=6)))
  }, simplify=FALSE))
  params <- list(thetas=initialThetas, alphas=initialAlphas,
                 deltas=initialDeltas, taus=initialTaus)
  chainMatrix <- matrix(0, nrow=length(unlist(params)), ncol=iterations)
  for (i in 1:iterations) {
    params <- iteration(params$thetas, params$alphas, params$deltas,
                        params$taus, responseMatrix)
    chainMatrix[ , i] <- unlist(params)
  }
  return(chainMatrix)
}

exampleMatrix <- ggum(responseMatrix, 3:5, 5000)
exampleMatrix2 <- ggum(responseMatrix, 3:5, 5000)
exampleMatrix3 <- ggum(responseMatrix, 3:5, 5000)

tracePlot <- function(trueValue, vectorList, paramName=''){
  vecLength <- unique(sapply(vectorList, length))
  if (length(vecLength) > 1) {stop('Vectors must be of same length', call.=F)}
  plot(1:vecLength, vectorList[[1]], type='l', col=rgb(0, 0, 0, alpha=0.75),
       main=paramName)
  sapply(2:length(vectorList), function(x){
    lineColor <- rgb(t(col2rgb(palette()[x])), alpha=191, maxColorValue=255)
    lines(1:vecLength, vectorList[[x]], col=lineColor)
    abline(h=mean(vectorList[[x]]), lty=x+2)
  })
  abline(h=trueValue, lty=2)
}

tracePlot(thetas[1], list(exampleMatrix[1, ], exampleMatrix2[1, ], exampleMatrix3[1, ]))
paramVec <- c(thetas, alphas, deltas, unlist(taus))
paramNames <- c(paste0('Theta', 1:10), paste0('Alpha', 1:3),
                paste0('Delta', 1:3), paste0('Tau 1-', 1:3),
                paste0('Tau 2-', 1:4), paste0('Tau 3-', 1:5))
invisible(sapply(1:28, function(x) tracePlot(paramVec[x], list(exampleMatrix[x, ],
                                                     exampleMatrix2[x, ],
                                                     exampleMatrix3[x, ]),
                                   paramNames[x])))


