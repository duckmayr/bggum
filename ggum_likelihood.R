# 
# We are going to do it like we did in cat, so it is exp(sum(log(P))) 
# rather than the product of the raw probabilities.  
# 
# It needs to take in as arguments:
#   
# 1) a matrix that has three columns:  
#       (1) the index for the respondent
#       (2) an index for the item 
#       (3) the observed answer.
# 2) A vector of thetas for the respondents
# 4) A list of lists containing the item parameters for the questions 
# (i think number of categories should just be one of the parameters that is passed in here).
# 
# For each row in the matrix, it needs to calculate the probability of that response given the 
# appropriate parameters (theta and item parameters).  The result of this should be vector.
# 
# Take the natural log of each element in the vector and sum them.
# 
# Take the exponent
# 
# Return


## probability<-function(k, numCats, theta, alpha, delta, tau)
probability<-function(k, numCats, theta, alpha, delta, tau){
  
  singleExp<-function(l, numCats, theta, alpha, delta, tau){
    return(exp(alpha*(l*(theta-delta) - sum(tau[1:(l+1)]))) +
             exp(alpha*((2*numCats-1-l)*(theta-delta)-sum(tau[1:(l+1)])))
    )
  }
  
  denom_vals = sapply(c(0:(numCats-1)), singleExp, numCats=numCats, theta=theta, alpha=alpha, delta=delta, tau=tau)
  numerator = denom_vals[k+1] 
  ## k takes on values {0,1,...,numCats-1}, but denom_vals indexes from {1,2,...,numCats}
  return (numerator / sum(denom_vals))
}

likelihood = function(response_matrix, theta_vec, item_params, mystery_alpha){
  probs = apply(response_matrix,1, function(row){
    ## grabbing info from the arguments supplied
    
    i = row[1]
    
    j = row[2]
    ans = row[3]
    
    theta = theta_vec[i]
    alpha = item_params$alpha[[j]]
    delta = item_params$delta[[j]]
    tau = item_params$tau[[j]]
    numCats = item_params$numCats[[j]]
    ## calling probability function
    row_prob = probability(ans, numCats, theta, alpha, delta, tau)
    
    return(row_prob)
  })
  probs[13] = probability(k = 3, numCats=5, theta=0, alpha=mystery_alpha, delta=-.5, tau= c(0,-1,-.7,-.5,-.2))
  ## log of probability vector
  log_probs = log(probs)
  ## exponent of the sum of the log probs
  return(exp(sum(log_probs)))
}

response_matrix = cbind(c(sapply(c(1:5),function(x){rep(x,5)})), ## respondents
                        rep(c(1:5), 5), ## items
                        sample(c(0:4), 25, replace=T)) ## responses

item_params = list("alpha"=c(runif(5,0,2)), "delta"=c(runif(5,0,1)), 
                   "tau"=c(lapply(c(1:5),function(x){return(c(0,sort(runif(4,-1,0))))})), "numCats"=c(rep(5,5)))
theta_vec = rnorm(5)

likelihood(response_matrix, theta_vec, item_params,mystery_alpha = 1)


