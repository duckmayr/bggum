#' Numerically stable product
#' 
#' \code{stableProd} calculates the product of a numeric vector
#' 
#' While base R provides the function \code{prod}, a more numerically stable
#' way to calculate the product of a vector is to take the exponential of the
#' sum of the natural logarithm of each element of the vector.
#' 
#' @param x A numeric vector
#' 
#' @return The product of the elements of x.
#' @export
stableProd <- function(x){
  return(exp(sum(log(x))))
}
