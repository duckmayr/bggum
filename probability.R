#' GGUM Probability Function
#' 
#' \code{probability} gives the probability a person will choose an item option.
#' 
#' Given parameters for the individual and the test item, \code{probabiliy}
#' returns the probability the individual will choose a particular response
#' option for the test item according to the GGUM.
#' 
#' @k A numeric vector; the index of the item option of interest
#' @K A numeric vector of length one; the number of options for the item
#' @theta A numeric vector of length one; the person's latent trait parameter
#' @alpha A numeric vector of length one; the item's discrimination parameter
#' @delta A numeric vector of length one; the item's location parameter
#' @tau  A numeric vector of length K; the vector threshold parameters for the
#'   item's options (where the first element of tau should be zero).
#' 
#' @return The probability that resopndent will choose response category k.
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