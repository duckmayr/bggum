#' GGUM Likelihood By Rows
#' 
#' \code{llRow} gives the (log) likelihood for one row of a response matrix
#' 
#' Given parameters for the individuals and test items, a response vector
#' recording an individual's response to each test item, \code{llRow} gives
#' the product of the probability corresponding to each entry in the relevant
#' vector.
#' 
#' @param theta A numeric vector of length one giving the individual's latent
#'   trait parameter
#' @param responseVector A numeric vector of length N (the number of
#'   respondents) giving the option chosen the individual for each item j
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
#' @param deltas A numeric vector of length n; each element of the vector is an
#'   item's location parameter
#' @param taus A list of numeric vectors; each list element j is a numeric
#'   vector of threshold parameters for item j's options (where the first
#'   element of the vector should be zero).
#' 
#' @return The (log) likelihood of the vector of interest.
#' @export
llRow <- function(theta, responseVector, alphas, deltas, taus){
      return(sum(log(sapply(1:length(responseVector), function(x){
         probability(responseVector[x], length(taus[[x]]), thetas, alphas[x],
                     deltas[x], taus[[x]])
      }))))
}

