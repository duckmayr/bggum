#' GGUM Likelihood By Columns
#' 
#' \code{llCol} gives the (log) likelihood for one column of a response matrix
#' 
#' Given parameters for the individuals and test items, a response vector
#' recording each individual's response to the test item, \code{llCol} gives
#' the product of the probability corresponding to each entry in the relevant
#' vector.
#' 
#' @param thetas A numeric vector of length N (the number of respondents); each
#'   each element of the vector is an individual's latent trait parameter
#' @param responseVector A numeric vector of length n (the number of items)
#'   giving the option chosen by individual i for the item
#' @param alpha A numeric vector of length one giving the item's discrimination
#'   parameter
#' @param delta A numeric vector of length one giving the item's location
#'   parameter
#' @param taus A numeric vector of length K (the number of options) giving the
#'   threshold parameter for each option; the first element should be zero.
#' 
#' @return The (log) likelihood of the vector of interest.
#' @export
llCol <- function(thetas, responseVector, alpha, delta, taus){
      return(sum(log(sapply(1:length(responseVector), function(x){
         probability(responseVector[x], length(taus), thetas[x], alpha, delta,
                     taus)
      }))))
}
