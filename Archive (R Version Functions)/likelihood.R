#' GGUM Likelihood
#' 
#' \code{likelihood} gives the likelihood of a response matrix or vector thereof
#' 
#' Given parameters for the individuals and test items, a response matrix
#' recording each individual's response to each test item, and as well as
#' (optionally) an indication of the vector of interest, \code{likelihood}
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
#' @param dimension A character vector of length one indicating whether the
#'   likelihood of a row or column should be given; should be one of 'row'
#'   or 'column'. If not given, the likelihood of the entire matrix will be
#'   calculated.
#' @param index A numeric vector of length one giving the row/column to be
#'   included in the calculation
#' 
#' @return The likelihood of the matrix or vector of interest.
#' @export
likelihood <- function(thetas, responseMatrix, alphas, deltas, taus,
                       dimension='', index=1){
  if ( dimension == 'row' ) {
    return(stableProd(sapply(1:ncol(responseMatrix), function(x){
      probability(responseMatrix[index, x], length(taus[[x]]), thetas[index],
                  alphas[x], deltas[x], taus[[x]])
    })))
  }
  if ( dimension == 'column' ) {
    return(stableProd(sapply(1:nrow(responseMatrix), function(x){
      probability(responseMatrix[x, index], length(taus[[index]]), thetas[x],
                  alphas[index], deltas[index], taus[[index]])
    })))
  }
  return(stableProd(probMatrix(thetas, responseMatrix, alphas, deltas, taus)))
}
