#' Calculate the mean silhouette score of a clustering.
#'
#' This function will compute the silhouette score for each cluster identified by `Seurat`'s Louvain modularity optimization community detection algorithm.
#'
#' @param seurat.obj The input object for which silhouette score will be computed. Defaults to NULL.
#' @param avg Should the average scores for each cluster be returned, or should a dataframe of every observation's cluster identity and score be returned? Defaults to TRUE.
#' @importFrom cluster silhouette
#' @export
#' @examples
#' ComputeSilhouetteScores(seurat.obj)

ComputeSilhouetteScores <- function(seurat.obj = NULL, avg = TRUE) {
  # check inputs
  if (is.null(seurat.obj)) stop ("Somehow, you didn't supply a Seurat object ...")
  # run function
  # prepare input matrix
  pca_df <- data.frame(Embeddings(seurat.obj, reduction = "pca"))
  pca_mat <- as.matrix(pca_df)
  # calculate cosine dissimilarity matrix
  cos_dist <- CosineDist(input = pca_mat)
  clust_list <- as.integer(seurat.obj$seurat_clusters) - 1
  res <- cluster::silhouette(dist = cos_dist, x = clust_list)
  # prepare results & return
  if (avg) {
    avg_widths <- summary(res)$clus.avg.widths
    avg_widths <- unlist(avg_widths)
    val <- avg_widths
  } else {
    val <- data.frame(Cluster = as.factor(res[, 1]),
                      Score = res[, 3])
  }
  return(val)
}
