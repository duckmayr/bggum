#' GGUM Likelihood By Rows
#' 
#' \code{likelihoodRow} gives the likelihood of the respondents (Rows)  of 
#' a response matrix or vector thereof
#' 
#' Given parameters for the individuals and test items, a response matrix
#' recording each individual's response to each test item, and as well as
#' (optionally) an indication of the vector of interest, \code{likelihoodRow}
#' gives the product of the probability corresponding to each entry in the
#' matrix or relevant vector.
#' 
#' @param thetas A numeric vector of length N (the number of respondents); each
#'   each element of the vector is an individual's latent trait parameter
#' @param responseMatrix A numeric matrix with N rows and n (the number of
#'   items) columns; each i,j element of the matrix gives the option chosen by
#'   individual i for item j
#' @param alphas A numeric vector of length n; each element of the vector is an
#'   item's discrimination parameter
#' @param deltas A numeric vector of length n; each element of the vector is an
#'   item's location parameter
#' @param taus A list of numeric vectors; each list element j is a numeric
#'   vector of threshold parameters for item j's options (where the first
#'   element of the vector should be zero).
#' @param index A numeric vector of length one giving the row to be
#'   included in the calculation
#' 
#' @return The likelihood of the matrix or vector of interest.
#' @export
likelihoodRow <- function(thetas, responseMatrix, alphas, deltas, taus, index=1){
      return(stableProd(sapply(1:ncol(responseMatrix), function(x){
         probability(responseMatrix[index, x], length(taus[[x]]), thetas[index],
                     alphas[x], deltas[x], taus[[x]])
      })))
}
