source('~/ggum/ggumR/R/probability.R') # change path as appropriate
source('~/ggum/ggumR/R/Beta_ab.R')

guessResponse <- function(theta, alpha, delta, tau){ # this is roughly the
  K <- length(tau) #                                 # process used in the
  probs <- probability(1:K, K, theta, alpha, delta, tau) # paper to make
  randomNumber <- runif(1) #                         # simulation data
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
thetas <- rnorm(400)
alphas <- rBeta_ab(10, params=alphaDistParams)
deltas <- rBeta_ab(10, params=deltaDistParams)
taus <- as.list(sapply(1:10, function(x){
  c(0, rBeta_ab(3, params=tauDistParams))
}, simplify=FALSE))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)

referenceFactor <- factor(c('strongly disagree', 'disagree',
                            'agree', 'strongly agree'),
                          levels=c('strongly disagree', 'disagree',
                                   'agree', 'strongly agree'))

quizData <- data.frame(Q1=rep('a', 400), Q2=rep('a', 400), Q3=rep('a', 400),
                       Q4=rep('a', 400), Q5=rep('a', 400), Q6=rep('a', 400),
                       Q7=rep('a', 400), Q8=rep('a', 400), Q9=rep('a', 400),
                       Q10=rep('a', 400))

for (i in 1:ncol(quizData)) {
  quizData[ , i] <- factor(levels(referenceFactor)[responseMatrix[ , i]],
                           levels=levels(referenceFactor))
}

all(sapply(quizData, as.numeric) == responseMatrix)

save(thetas, alphas, deltas, taus, quizData, responseMatrix,
     file='~/ggum/exampleData.RData')

## Make a larger data set:

thetas <- rnorm(1000)
alphas <- rBeta_ab(20, params=alphaDistParams)
deltas <- rBeta_ab(20, params=deltaDistParams)
taus <- as.list(sapply(1:20, function(x){
  c(0, rBeta_ab(3, params=tauDistParams))
}, simplify=FALSE))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)
apply(responseMatrix, 2, function(x) length(unique(x)))
save(thetas, alphas, deltas, taus, responseMatrix, file='lgExData.RData')

## While the last simulated dataset eliminated the unanimous response problem,
## we still had some questions where not all options were selected. In order
## to estimate all of the true parameters, we need a dataset where all options
## are selected, so I made sure every option was picked for every question:

thetas <- rnorm(1000)
alphas <- rBeta_ab(20, params=alphaDistParams)
deltas <- rBeta_ab(20, params=deltaDistParams)
taus <- as.list(sapply(1:20, function(x){
  c(0, rBeta_ab(3, params=tauDistParams))
}, simplify=FALSE))
responseMatrix <- mkResponseMatrix(thetas, alphas, deltas, taus)
apply(responseMatrix, 2, function(x) length(unique(x)))
# Only questions 1, 4, 12, and 15 have not all options selected.
# So, we make sure all options are selected:
q1 <- responseMatrix[ , 1]
length(unique(q1))
while ( length(unique(q1)) != 4 ) {
  q1 <- sapply(1:1000, function(x){
    guessResponse(thetas[x], alphas[1], deltas[1], taus[[1]])
  })
}
length(unique(q1))
q4 <- responseMatrix[ , 4]
length(unique(q4))
while ( length(unique(q4)) != 4 ) {
  q4 <- sapply(1:1000, function(x){
    guessResponse(thetas[x], alphas[4], deltas[4], taus[[4]])
  })
}
length(unique(q4))
q12 <- responseMatrix[ , 12]
length(unique(q12))
while ( length(unique(q12)) != 4 ) {
  q12 <- sapply(1:1000, function(x){
    guessResponse(thetas[x], alphas[12], deltas[12], taus[[12]])
  })
}
# deltas[12] was -3.8... impossible to get all responses
# (confirmed via plotting item response function)
# So, we have to replace this delta to be able to get all responses
deltas[12] <- rBeta_ab(1, params=deltaDistParams)
while ( length(unique(q12)) != 4 ) {
  q12 <- sapply(1:1000, function(x){
    guessResponse(thetas[x], alphas[12], deltas[12], taus[[12]])
  })
}
length(unique(q12))
q15 <- responseMatrix[ , 15]
length(unique(q15))
while ( length(unique(q15)) != 4 ) {
  q15 <- sapply(1:1000, function(x){
    guessResponse(thetas[x], alphas[15], deltas[15], taus[[15]])
  })
}
length(unique(q15))
responseMatrix[ , 1] <- q1
responseMatrix[ , 4] <- q4
responseMatrix[ , 12] <- q12
responseMatrix[ , 15] <- q15
apply(responseMatrix, 2, function(x) length(unique(x))) # all have all answers
save(thetas, alphas, deltas, taus, responseMatrix, file='newExData.RData')

