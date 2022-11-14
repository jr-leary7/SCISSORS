#' Calculate the silhouette score of a clustering.
#'
#' @name ComputeSilhouetteScores
#' @author Jack Leary
#' @description This function will compute the silhouette score for each cluster identified by `Seurat`'s Louvain modularity optimization community detection algorithm.
#' @importFrom Seurat Embeddings
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @param seurat.obj The input object for which silhouette score will be computed. Defaults to NULL.
#' @param dist.metric Which distance metric should be used? Defaults to "cosine", but any of the metrics used by \code{\link[stats]{dist}} will work.
#' @param avg Should the average scores for each cluster be returned, or should a dataframe of every observation's cluster identity and score be returned? Defaults to TRUE.
#' @return If \code{avg = TRUE}, returns the average silhouette score per cluster, else returns a cell-level dataframe of the cluster identities & silhouette scores.
#' @seealso \code{\link{CosineDist}}.
#' @export
#' @examples
#' \dontrun{
#' ComputeSilhouetteScores(seurat.obj)
#' ComputeSilhouetteScores(seurat.obj,
#'                         dist.metric = "euclidean",
#'                         avg = FALSE)
#' }

ComputeSilhouetteScores <- function(seurat.obj = NULL,
                                    dist.metric = "cosine",
                                    avg = TRUE) {
  # check inputs
  if (is.null(seurat.obj)) { stop("You didn't supply a Seurat object to ComputeSilhouetteScores().") }
  # run function
  # prepare input matrix
  pca_df <- data.frame(Seurat::Embeddings(seurat.obj, reduction = "pca"))
  pca_mat <- as.matrix(pca_df)
  # calculate distance matrix -- default is cosine dissimilarity
  if (dist.metric == "cosine") {
    pc_dists <- CosineDist(input.mat = pca_mat)
  } else {
    pc_dists <- stats::dist(x = pca_mat, method = dist.metric)
  }
  clust_list <- as.integer(seurat.obj$seurat_clusters) - 1
  res <- cluster::silhouette(dist = pc_dists, x = clust_list)
  # prepare results & return
  if (avg) {
    avg_widths <- summary(res)$clus.avg.widths
    avg_widths <- unlist(avg_widths)
    val <- avg_widths
  } else {
    val <- data.frame(Cluster = as.factor(res[, 1]), Score = res[, 3])
  }
  return(val)
}
