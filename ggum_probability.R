## Matt -- you had sum(tau[1:l+1])
## --- needed to be sum(tau[1:(l+1)])


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

probability(0, 5, 0, 2, 1, c(0,1,2, 3, 4))

### replicating figure 1
testFun<-function(x, k){
  return(probability(k, numCats = 4, theta = x, alpha = 2, delta = 0, tau = c(0,-1, -.7, -.4))) # (tau_0 = 0 by definition)
}
testFun(0, 3)

dev.off()
par(mfrow=c(1,1))
theta_vec<-seq(-3, 3, by=.05)
cols = c('red', 'blue', 'green', 'purple')
for(k in 0:3){
  if(k==0){
  plot(theta_vec, sapply(theta_vec, function(x){testFun(x,k)}), 
       xlab="theta", ylab = paste0("P(k|theta)"), ylim = c(0,1),
       cex=.5,col=cols[k+1])}
  else{
  lines(theta_vec, sapply(theta_vec, function(x){testFun(x,k)}), 
       xlab="", ylab ='', ylim = c(0,1),
       cex=.5,col=cols[k+1])
  }
  par(new=T)
} 
legend('topright', fill = cols, legend=paste("k =",c(0:3)),cex=.6)


?plot
