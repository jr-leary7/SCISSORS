#' A function to iteratively recalculate cluster assignments based on highly variable genes
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster.
#' @import Seurat
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param recluster.res The amount of times to recluster each cluster and subcluster. Defaults to 1.
#' @param var.method How should highly variable genes be identified? One of "sctransform" or "vst". Defaults to "sctransform".
#' @param n.variable.genes How many variable genes should be detected in each subcluster? Defaults to 3000
#' @param use.SingleR.labels Should cell types from `SingleR` be used to recluster cells? Defaults to FALSE.
#' @export
#' @examples
#' ReclusterCells(seurat.object)
#' ReclusterCells(seurat.object, recluster.res = 2, use.SingleR.labels = TRUE)

ReclusterCells <- function(seurat.object = NULL,
                           recluster.res = 1,
                           var.method = "sctransform",
                           use.SingleR.labels = FALSE) {
  # test input
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object as input!") }
  # run on clusters assigned by Seurat
  if (!use.SingleR.labels) {
    reclust_list <- list()
    for (clust in seq(seurat.object$seurat_clusters)) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      if (var.method == "sctransform") {
        temp_obj <- SCTransform(temp_obj,
                                variable.features.n = n.variable.genes,
                                seed.use = 629,
                                verbose = FALSE)
      } else if (var.method == "vst") {
        temp_obj <- FindVariableFeatures(temp_obj,
                                         selection.method = "vst",
                                         nfeatures = n.variable.genes)
      } else { stop("Choose a viable highly variable gene selection method") }

      # run PCA
      temp_obj <- RunPCA(temp_obj,
                         npcs = 1:20,
                         verbose = FALSE,
                         seed.use = 629)
      # find clusters
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                dims = 1:20)
      temp_obj <- FindClusters(temp_obj,
                               resolution = .3,
                               random.seed = 629)

      reclust_list[[clust]] <- NULL  # change this to be the correct output later
    }
  }
}
