#' A function to compute the cosine similarity of two vectors.
#'
#' This function takes as input two vectors of length *n*, *x* and *y*, for which the cosine similarity is computed.
#' @param x Input vector. Default is NULL.
#' @param y Input vector. Default is NULL
#' @examples
#' CosSimilarity(x, y)

CosSimilarity <- function(x = NULL, y = NULL) {
  # test inputs
  if (is.null(x) | is.null(y)) { stop("You must provide two vectors as input!") }
  if (length(x) != length(y)) { stop("Input vectors must have the same number of elements!") }

  # calculate cosine similarity
  cos_sim <- sum((x * y)) / ((sqrt(sum(x^2))) * (sqrt(sum(y^2))))
  return(cos_sim)
}
