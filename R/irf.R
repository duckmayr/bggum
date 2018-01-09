#' Item Response Function
#' 
#' Plots a response function given alpha, delta, and tau parameters.
#' 
#' @param a A numeric vector of length one, the alpha parameter
#' @param d A numeric vector of length one, the delta parameter
#' @param t A numeric vector of length K (the number of options),
#'   the vector of tau parameters for each option --
#'   note the first element of the vector should be zero
#' @param from A numeric vector of length one,
#'   the lowest theta value to estimate response probabilities for
#' @param to A numeric vector of length one,
#'   the highest theta value to estimate response probabilities for
#' @param by A numeric vector of length one giving the spacing between
#'   theta values
#' @param sub A subtitle for the resulting plot
#' 
#' @export
irf <- function(a, d, t, from=-3, to=3, by=0.01, sub=''){
  th <- seq(from=from, to=to, by=by)
  K <- length(t)
  plot(th, sapply(th, function(x) ggumProbability(1, x, a, d, t)),
       type='l', xlab=expression(theta), ylab=expression(P[ij](k)),
       main=paste('Item Response Function', sub),
       xlim=c(from, to), ylim=c(0, 1))
  for ( i in 2:K ) {
    lines(th, sapply(th, function(x) ggumProbability(i, x, a, d, t)), lty=i)
  }
}
