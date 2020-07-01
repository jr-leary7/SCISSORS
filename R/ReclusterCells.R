#' A function to iteratively recalculate cluster assignments based on highly variable genes.
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using the variance-stabilizing transformation implemented in `Seurat`.
#' @import Seurat
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param n.variable.genes How many variable genes should be detected in each subcluster? Defaults to 3000
#' @param use.SingleR.labels Should cell types from `SingleR` be used to recluster cells? Defaults to FALSE.
#' @param which.clust How should we decide which clusters to recluster? Defaults to "auto", which tests homogeneity of clusters, and reclusters the heterogenous ones. Can also be a user-defined list of clusters, which overrides the cluster homogeneity test.
#' @param random.seed The random seed used to control stochasticity. Defaults to 629.
#' @export
#' @examples
#' ReclusterCells(seurat.object)
#' ReclusterCells(seurat.object, n.variable.genes = 4000, use.SingleR.labels = TRUE)

ReclusterCells <- function(seurat.object = NULL,
                           use.SingleR.labels = FALSE,
                           n.variable.genes = 3000,
                           which.clust = "auto") {
  # test input
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object as input!") }
  reclust_list <- list()
  if (! use.SingleR.labels & which.clust == "auto") {
    print("Identifying most variable genes and re-clustering each cluster")
    unique_clusts <- sort(as.integer(unique(seurat.object$seurat_clusters)) - 1)
    for (clust in unique_clusts) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      temp_obj <- FindVariableFeatures(temp_obj,
                                       assay = "SCT",
                                       selection.method = "vst",
                                       nfeatures = n.variable.genes,
                                       verbose = FALSE)
      temp_obj <- RunPCA(temp_obj,
                         npcs = 30,
                         seed.use = 629,
                         verbose = FALSE)
      # use k ~ sqrt(N) general rule
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                k.param = round(sqrt(ncol(temp_obj))))
      # cluster several times, compute silhouette scores, compare & choose best re-clustering
      sil_scores <- c()
      res_vals <- c(.075, .2, .35)
      for (res in seq(res_vals)) {
        temp_obj <- FindClusters(temp_obj,
                                 resolution = .1,
                                 algorithm = 1,
                                 random.seed = random.seed,
                                 verbose = FALSE)
        sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
        mean_sil <- mean(sil_res)
        sil_scores[res] <- mean_sil
      }
      names(sil_scores) <- c(".075", ".2", ".35")
      correct_res <- as.integer(names(sil_scores[sil_scores == max(sil_scores)]))
      temp_obj <- FindClusters(temp_obj,
                               resolution = correct_res,
                               algorithm = 1,
                               random.seed = random.seed)

    }
    names(reclust_list) <- unique_clusts
  }

  if (! use.SingleR.labels & which.clust == "simple") {
    print("Identifying most variable genes and re-clustering each cluster")
    unique_clusts <- sort(as.integer(unique(seurat.object$seurat_clusters)) - 1)
    for (clust in unique_clusts) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      temp_obj <- FindVariableFeatures(temp_obj,
                                       assay = "SCT",
                                       selection.method = "vst",
                                       nfeatures = n.variable.genes,
                                       verbose = FALSE)
      temp_obj <- RunPCA(temp_obj,
                         npcs = 30,
                         seed.use = 629,
                         verbose = FALSE)
      # use k ~ sqrt(N) general rule
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                k.param = round(sqrt(ncol(temp_obj))))
      temp_obj <- FindClusters(temp_obj,
                               resolution = .1,
                               algorithm = 1,
                               random.seed = random.seed)
      reclust_list[[clust]] <- temp_obj
    }
    names(reclust_list) <- unique_clusts
  }
  # merge re-clustered results back in to original Seurat object
  max_clust_id <- max(unique_clusts)
  for (clust in unique_clusts) {
    temp_obj <- reclust_list[[clust]]
    if (length(unique(temp_obj$seurat_clusters)) > 1) {

    }
  }

  return(NULL)  # change later
}
