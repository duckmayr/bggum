#' probability
#' 
#' @k Option k of item j
#' @numCats The number of response options for the item: k = {0,1,...,numCats-1}
#' @theta The postion of the respondent on the latent trait scale
#' @alpha The discrimination parameter of item j
#' @delta The location parameter of item j
#' @tau  A vector of K-1 threshold parameters: {tau_j1, tau_j2, ..., tau_j(k-1)}
#' 
#' @return The probability that resopndent will choose response category k.
probability<-function(k, numCats, theta, alpha, delta, tau){
                
      denom_vals = sapply(c(0:(numCats-1)), function(l){
          print(l)
          return(exp(alpha*(l*( (theta-delta))-sum(tau))) +
              exp(alpha*( (2*numCats-1-l)*(theta-delta)-sum(tau) )))
      })
      numerator = denom_vals[k+1] 
      print(denom_vals)
      print(numerator)
      ## k takes on values {0,1,...,numCats-1}, but denom_vals indexes from {1,2,...,numCats}
      return (numerator / sum(denom_vals))
}

probability(2, 5, 5, 2, 1, c(-1, -2, 3, 4))

testFun<-function(x, k){
  probability(k, 5, x, 2, 1, c(-1, -2, 3, 4))
}

x.axis<-seq(-10, 10, by=.1)
plot(sapply(x.axis, testFun, k=4), sapply(x.axis, testFun, k=2))

