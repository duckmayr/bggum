#' Error-free (potential) extraction
#' 
#' \code{safeExtract} extracts an element of a vector if and only if it exists.
#' 
#' @param vector The vector from which to extract an element
#' @param index The index of the element to extract
#' 
#' @return If the vector has at least index elements, the index element of
#'   vector is returned; otherwise, nothing is returned.
safeExtract <- function(vector, index){
  if (length(vector) >= index) {
    return(unlist(vector[index]))
  }
}