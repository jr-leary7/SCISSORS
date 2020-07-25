#' Compute the cosine distance matrix.
#'
#' This function takes a matrix as input, and computes the cosine distance (1 - cosine similarity) between the observations.
#' @param input The input matrix. Defaults to NULL.
#' @export
#' @examples
#' CosineDist(input = pca_matrix)

CosineDist <- function(input = NULL) {
  # check inputs -- although as this is a helper function, it should never be called incorrectly
  if (is.null(input)) { stop("You forgot to provide an input matrix") }
  # run function
  dist_mat <- as.dist(1 - input %*% t(input) / (sqrt(rowSums(input^2) %*% t(rowSums(input^2)))))
  return(dist_mat)
}
