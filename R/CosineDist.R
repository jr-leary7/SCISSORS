#' Compute the cosine distance matrix.
#'
#' @name CosineDist
#' @author Jack Leary
#' @description This function takes a matrix as input, and computes the cosine distance (1 - cosine similarity) between the observations.
#' @importFrom stats as.dist
#' @importFrom coop tcosine
#' @param input.mat The input matrix. Defaults to NULL.
#' @export
#' @examples
#' \dontrun{
#' CosineDist(input.mat = pca_matrix)
#' }

CosineDist <- function(input.mat = NULL) {
  # check inputs
  if (is.null(input.mat)) { stop("You must provide a matrix to CosineDist().") }
  # compute cosine distance
  dist_mat <- stats::as.dist(1 - coop::tcosine(x = input.mat))
  return(dist_mat)
}
