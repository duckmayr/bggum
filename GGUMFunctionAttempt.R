#' probability
#' 
#' @answer The response chosen by the respondent
#' @numCats The number of response options for the item
#' @theta The postion of the respondent on the latent trai
#' @alpha The discrimination parameter
#' @delta The location parameter
#' @tau  A vector of K-1 threshold parameter
#' 
#' @return The probability that resopndent as position theta will choose category answer for this item.
probability<-function(answer, numCats, theta, alpha, delta, tau){
                      return(exp(alpha*(answer*(theta-delta))-sum(tau))+exp(alpha*(2*numCats-1-answer)*(theta-delta)-sum(tau)))
}

probability(2, 5, 50, 2, 1, c(-1, -2, 3, 4))

testFun<-function(x, k){
  probability(k, 5, x, 2, 1, c(-1, -2, 3, 4))
}

x.axis<-seq(-10, 10, by=.1)
plot(sapply(x.axis, testFun, k=4), sapply(x.axis, testFun, k=2))

