#' A function to compute the cosine distance matrix.
#'
#' This function takes a matrix as input, and computes the cosine distance (1 - cosine similarity) between the observations.
#' @param input The input matrix. Defaults to NULL.
#' @examples
#' CosineDist(input = pca_matrix)

CosineDist <- function(input = NULL) {
  if (is.null(input)) { stop("You forgot to provide an input matrix") }
  dist_mat <- as.dist(1 - input %*% t(input) / (sqrt(rowSums(input^2) %*% t(rowSums(input^2)))))
  return(dist_mat)
}
