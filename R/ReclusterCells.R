#' A function to iteratively recalculate cluster assignments based on highly variable genes.
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using the variance-stabilizing transformation implemented in `Seurat`.
#' @import Seurat
#' @importFrom ggplot2 labs
#' @importFrom plyr mapvalues
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param n.variable.genes How many variable genes should be detected in each subcluster? Defaults to 3000
#' @param which.clust How should we decide which clusters to recluster? Defaults to "auto", which judges the quality of a re-clustering using silhouette scoring. Can also be a user-defined vector of clusters to investigate.
#' @param resolution.vals (Optional) A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.05, .1, .15, .2, .35).
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' ReclusterCells(seurat.object)
#' ReclusterCells(seurat.object, n.variable.genes = 4000, use.SingleR.labels = TRUE)

ReclusterCells <- function(seurat.object = NULL,
                           n.variable.genes = 4000,
                           which.clust = NULL,
                           resolution.vals = c(.05, .1, .2, 4)) {
  # test input
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object as input!") }

  # auto-recluster cells on a per-cluster basis using silhouette scores to evaluate result
  reclust_list <- list()
  res_df <- data.frame(CellID = rownames(seurat.object@meta.data),
                       ClustID = rep(NA, nrow(seurat.object@meta.data)))
  if (which.clust == "auto") {
    print(sprintf("Identifying %s most variable genes and re-clustering each cluster", n.variable.genes))
    unique_clusts <- sort(as.integer(unique(seurat.object$seurat_clusters)) - 1)
    for (clust in unique_clusts) {
      # reprocess each cluster as its own object
      temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      temp_obj <- SCTransform(temp_obj,
                              vars.to.regress = "percent_MT",
                              seed.use = random.seed,
                              variable.features.n = n.variable.genes,
                              verbose = FALSE)
      temp_obj <- RunPCA(temp_obj,
                         npcs = 15,
                         seed.use = 629,
                         verbose = FALSE,
                         features = VariableFeatures(temp_obj))
      temp_obj <- RunTSNE(temp_obj,
                          reduction = "pca",
                          dim.embed = 2,
                          seed.use = random.seed)
      # use k ~ sqrt(N) general rule to recreate the SNN graph
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                k.param = round(sqrt(ncol(temp_obj))))
      # cluster several times, compute silhouette scores, compare & choose best re-clustering
      sil_scores <- c()
      for (res in seq(resolution.vals)) {
        temp_obj <- FindClusters(temp_obj,
                                 resolution = resolution.vals[res],
                                 algorithm = 1,
                                 random.seed = random.seed,
                                 verbose = FALSE)
        if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
          sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
          mean_sil <- mean(sil_res)
          sil_scores[res] <- mean_sil
        } else {
          # neutral placeholder value for the case when the number of identified clusters is 1
          sil_scores[res] <- 0
        }
      }
      names(sil_scores) <- as.character(resolution.vals)
      # save correct results using resolution that gives highest mean silhouette score
      if (max(sil_scores) > .25) {
        correct_res <- as.numeric(names(sil_scores[sil_scores == max(sil_scores)]))
        print(sprintf("Clustering cells using resolution = %s, which achieved silhouette score: %s",
                      correct_res,
                      round(max(sil_scores), 4)))
        temp_obj <- FindClusters(temp_obj,
                                 resolution = correct_res,
                                 algorithm = 1,
                                 random.seed = random.seed)
      } else {
        # replace new object w/ original one, & find variable features again just in case
        print(sprintf("Did not find suffcient evidence of subclusters, as the max silhouette score was: %s", round(max(sil_scores), 4)))
        temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      }
      reclust_list[[clust + 1]] <- temp_obj
    }
    names(reclust_list) <- unique_clusts
    # reset cluster IDs so as to have unique values
    max_clust_id <- 0
    for (clust in unique_clusts) {
      temp_obj <- reclust_list[[clust + 1]]
      clust_ids <- as.integer(temp_obj$seurat_clusters) - 1
      res_df[res_df$CellID %in% colnames(temp_obj)]$ClustID <- clust_ids
      max_clust_id <- max(res_df$ClustID[!is.na(res_df$ClustID)])
    }

  } else if (!is.null(which.clust)) {
    print(sprintf("Identifying subclusters using %s variable genes in clusters: %s"), n.variable.genes, which.clust)
    for (clust in seq(which.clust)) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[clust])
      temp_obj <- FindVariableFeatures(temp_obj,
                                       assay = "SCT",
                                       selection.method = "vst",
                                       verbose = FALSE,
                                       nfeatures = n.variable.genes)
      temp_obj <- RunPCA(temp_obj,
                         features = VariableFeatures(temp_obj),
                         npcs = 15,
                         verbose = FALSE,
                         seed.use = random.seed)
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                dims = 1:15,
                                k.param = round(sqrt(ncol(temp_obj))),
                                verbose = FALSE)
      sil_scores <- c()
      for (res in seq(resolution.vals)) {
        temp_obj <- FindClusters(temp_obj,
                                 resolution = resolution.vals[res],
                                 algorithm = 1,
                                 random.seed = random.seed,
                                 verbose = FALSE)
        print(DimPlot(temp_obj))
        if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
          sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
          mean_sil <- mean(sil_res)
          sil_scores[res] <- mean_sil
        } else {
          # neutral placeholder value for the case when the number of identified clusters is 1
          sil_scores[res] <- 0
        }
      }
      # save best results
      names(sil_scores) <- as.character(resolution.vals)
      if (max(sil_scores) > .2) {
        correct_res <- as.numeric(names(sil_scores[sil_scores == max(sil_scores)]))
        temp_obj <- FindClusters(temp_obj,
                                 resolution = correct_res,
                                 algorithm = 1,
                                 random.seed = random.seed)
      } else {
        # replace new object w/ original one, & find variable features again just in case
        temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
        temp_obj <- FindVariableFeatures(temp_obj,
                                         assay = "SCT",
                                         selection.method = "vst",
                                         nfeatures = n.variable.genes,
                                         verbose = FALSE)
      }
    }
    # merge results and return
  }
}




