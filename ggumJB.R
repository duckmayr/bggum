
probability <- function(k, K, theta, alpha, delta, tau){
  # Equation 1 in de la Torre 2006 is a quotient where the numerator expression
  # is summed for the denominator expression. The numerator expression can be
  # further broken down into two smaller expressions, each of which uses a
  # common espression.
  smExpr <- cumsum(tau) # expression common to both numerator expressions
  lgExpr1 <- exp(alpha * (0:(K-1) * (theta - delta) - smExpr))
  lgExpr2 <- exp(alpha * ((2*K - (0:(K-1)) - 1) * (theta - delta) - smExpr))
  numerator <- lgExpr1 + lgExpr2
  return(numerator[k]/sum(numerator)) # we only want the kth element
}


# So let's try out the function and replicate Figure 1 from de la Torre 2006
probability(1, 5, 0, 2, 1, c(0, 1, 2, 3, 4))
a <- 2
d <- 0
t <- c(0, -1, -0.7, -0.4)
th <- seq(from=-3, to=3, by=0.1)
probs <- t(sapply(th, function(x){
  sapply(1:length(t), probability, length(t), x, a, d, t)
}))
plot(NULL, ylim=c(0, 1), xlim=c(-3, 3), ylab='Prob. of Positive Response',
     xlab=expression(theta), main='GGUM Option Response Functions')
plotColors <- c('grey40', 'grey50', 'grey20', 'grey30')
plotSymbols <- c(18, 20, 17, 15)
indices <- c(TRUE, FALSE)
matlines(th, probs, col=plotColors, lty=1)
matpoints(th[indices], probs[indices, ], col=plotColors, pch=plotSymbols)
legend('right', legend=c('Opt. 1', 'Opt. 2', 'Opt. 3', 'Opt. 4'),
       pch=plotSymbols, col=plotColors)


probMatrix <- function(thetas, responseMatrix, alphas, deltas, taus){
  # Given vectors (or for tau, a list) of parameters, we can calculate the
  # probability of each choice in the reponse matrix
  return(t(sapply(1:length(thetas), function(x){ # for each person
    sapply(1:length(alphas), function(y){ # for each item
      probability(responseMatrix[x, y], length(taus[[y]]), thetas[x],
                  alphas[y],  deltas[y], taus[[y]])})})))
}

stableProd <- function(x){
  # More numerically stable than prod() for likelhood of the response matrix
  return(exp(sum(log(x))))
}

likelihood <- function(thetas, responseMatrix, alphas, deltas, taus){
  # We just get the product of all the probabilities for the response matrix
  return(stableProd(probMatrix(thetas, responseMatrix, alphas, deltas, taus)))
}

thetas <- rnorm(100)
responses <- matrix(c(sample(1:3, 100, replace=TRUE),
                      sample(1:4, 100, replace=TRUE),
                      sample(1:5, 100, replace=TRUE)), nrow=100)
alphas <- rnorm(3)
deltas <- rnorm(3)
taus <- list(c(0, rnorm(2)),c(0, rnorm(3)), c(0, rnorm(4)))

likelihood(thetas, responses, alphas, deltas, taus)

## To make some fake data using the most likely response given parameters:
likelyResponse <- function(theta, alpha, delta, tau){
  K <- length(tau)
  return(which.max(sapply(1:K, probability, K, theta, alpha, delta, tau)))
}

likelyResponseVector <- function(theta, alpha, delta, tau){
  return(sapply(1:length(tau), function(x){
    likelyResponse(theta, alpha[x], delta[x], tau[[x]])}))
}

plotLcurve <- function(trueTheta, response, alpha, delta, tau){
  xVals <- seq(from=-3, to=3, by=0.1)
  yVals <- sapply(xVals, likelihood, response, alpha, delta, tau)
  plot(x=xVals, y=yVals, type='l', xlab=expression(theta), ylab='Likelihood',
       main=expression(paste('Likelihood of Response at Values of ', theta)))
  abline(v=trueTheta, lty=2)
}

trueThetas <- c(-3, -2, -1, 0, 1, 2, 3)
alpha <- c(1, 2, 3)
delta <- c(3, 2, 1)
tau <- list(c(0, -1, -2, -3), c(0, -3, -2, -1), c(0, -2, -1, -3))
invisible(sapply(trueThetas, function(x){
  plotLcurve(x, matrix(likelyResponseVector(x, alpha, delta, tau), nrow=1),
             alpha, delta, tau)
}))

#######################

rBeta_ab <- function(n, shape1=2, shape2=3, a = 0, b = 1,
                     params = list(shape1, shape2, a, b),...){
  if(!missing(params)){
    shape1 <- params$shape1
    shape2 <- params$shape2
    a <- params$a
    b <- params$b
  }
  X <- rbeta(n,shape1,shape2)
  out <- (b-a)*X + a
  return(out)
}

dBeta_ab <-function(x, shape1=2, shape2=3, a = 0, b=1,
                    params = list(shape1, shape2, a, b),...){
  if(!missing(params)){
    shape1 <- params$shape1
    shape2 <- params$shape2
    a <- params$a
    b <- params$b
  }
  out <- (x>=a & x<=b) * dbeta((x-a)/(b-a),shape1,shape2)/(b-a)
  return(out)
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
responseMat
likelihood(thetas, responseMat, alphas, deltas, taus)

getPrior <- function(param, value){
  alphaDistParams <- list(shape1=1.5, shape2=1.5, a=.25, b=4)
  deltaDistParams <- list(shape1=2, shape2=2, a=-5, b=5)
  tauDistParams <- list(shape1=2, shape2=2, a=-6, b=6)
  return(switch(param, 'alpha'=dBeta_ab(value, params=alphaDistParams),
                'delta'=dBeta_ab(value, params=deltaDistParams),
                'tau'=dBeta_ab(value, params=tauDistParams),
                'theta'=dnorm(value)))
}

ratio <- function(param, person, item, currentValue, proposedValue, alphas,
                  deltas, taus, thetas, responseMatrix, paramType){
  probMat <- probMatrix(thetas, responseMatrix, alphas, deltas, taus)
  currentLikelihood <- stableProd(probMat)
  # paramType <- switch(param+1, 'theta', 'alpha', 'delta', 'tau')
  if (paramType == 'theta') {
    probMat[person, ] <- sapply(1:ncol(probMat), function(x){
      probability(responseMatrix[person, x], length(taus[[x]]), proposedValue,
                  alphas[x], deltas[x], taus[[x]])
    })
  } else {
    params <- unlist(c(alphas[item], deltas[item], taus[[item]]))
    params[param] <- proposedValue
    probMat[ , item] <- sapply(1:nrow(probMat), function(x){
      probability(responseMatrix[x, item], length(taus[[item]]), thetas[x],
                  params[1], params[2], params[3:length(params)])
    })
  }
  proposedLikelihood <- stableProd(probMat)
  return((proposedLikelihood*getPrior(paramType, proposedValue))
         /(currentLikelihood*getPrior(paramType, currentValue)))
}

ratio(0, 1, 1:length(alphas), thetas[1], -1, alphas, deltas, taus, thetas, responseMat, 'theta')
ratio(0, 1, 1:length(alphas), thetas[1], 0, alphas, deltas, taus, thetas, responseMat, 'theta')
ratio(0, 1, 1:length(alphas), thetas[1], 1, alphas, deltas, taus, thetas, responseMat, 'theta')
ratio(1, 1:length(thetas), 1, alphas[1], 1, alphas, deltas, taus, thetas, responseMat, 'alpha')
ratio(2, 1:length(thetas), 1, deltas[1], 1, alphas, deltas, taus, thetas, responseMat, 'delta')
ratio(4, 1:length(thetas), 1, taus[[1]][2], 1, alphas, deltas, taus, thetas, responseMat, 'tau')

proposer <- function(currentValue, paramType, param, paramVector, item){
  if (paramType != 'tau') {
    return(rnorm(1, currentValue, sd(paramVector)))
  }
  paramVector <- sapply(paramVector, '[[', param-2)
  return(rnorm(1, currentValue, sd(paramVector)))
}

iteration <- function(param, person, item, currentValue, alphas, deltas, taus,
                      thetas, responseMatrix, paramType){
  proposedValue <- proposer(currentValue, paramType, param,
                            switch(paramType, 'theta'=thetas, 'alpha'=alphas,
                                   'delta'=deltas, 'tau'=taus), item)
  reSetProb <- ratio(param, person, item, currentValue, proposedValue, alphas,
                     deltas, taus, thetas, responseMatrix, paramType)
  if (reSetProb > 1) {reSetProb <- 1}
  reSet <- sample(c(TRUE, FALSE), 1, prob=c(reSetProb, 1-reSetProb))
  return(ifelse(reSet, proposedValue, currentValue))
}

iteration(0, 1, 1:length(alphas), thetas[1], alphas, deltas, taus, thetas, responseMat, 'theta')
iteration(0, 1, 1:length(alphas), thetas[1], alphas, deltas, taus, thetas, responseMat, 'theta')
iteration(0, 1, 1:length(alphas), thetas[1], alphas, deltas, taus, thetas, responseMat, 'theta')
iteration(1, 1:length(thetas), 1, alphas[1], alphas, deltas, taus, thetas, responseMat, 'alpha')
iteration(2, 1:length(thetas), 1, deltas[1], alphas, deltas, taus, thetas, responseMat, 'delta')
iteration(4, 1:length(thetas), 1, taus[[1]][2], alphas, deltas, taus, thetas, responseMat, 'tau')
