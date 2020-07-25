#' Identify subpopulations in single cell clusters.
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using `sctransform` as implemented in `Seurat`. The function returns a list of `Seurat` objects, one for each cluster the user wants to investigate.
#' @import Seurat
#' @importFrom ggplot2 labs
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param n.variable.genes How many variable genes should be detected in each subcluster? Defaults to 4000.
#' @param n.PC How many PCs should be used as input to non-linear to non-linear dimension reduction and clustering algorithms. Defaults to 10.
#' @param redo.tSNE (Optional) Should a cluster-specific t-SNE embedding be generated? Sometimes subpopulations appear mixed together on the original t-SNE coordinates, but separate clearly when re-embedded. Defaults to TRUE.
#' @param which.clust Which clusters should undergo subpopulation detection analysis? If a user-defined list is not provided, all clusters will be re-clustered; this is NOT recommended as some clusters will not contain subpopulations.
#' biologically non-relevant clusters will be "discovered," and the FP rate will increase. Defaults to "auto".
#' @param resolution.vals (Optional) A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.05, .1, .15, .2, .35).
#' @param k.val (Optional) The parameter *k* to be used when creating the shared nearest-neighbor graph. Defaults to *k* ~ sqrt(*n*).
#' @param do.plot (Optional) Should t-SNE plots of the various reclusterings be plotted for visual inspection by the user? Defaults to FALSE.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' ReclusterCells(seurat.object)
#' ReclusterCells(seurat.object, resolution.vals = c(0.1, 0.2, 0.3), random.seed = 100)
#' ReclusterCells(seurat.object, n.variable.genes = 3000, which.clust = list(0, 3, 5), do.plot = TRUE)

ReclusterCells <- function(seurat.object = NULL,
                           n.variable.genes = 4000,
                           n.PC = 10,
                           redo.tSNE = TRUE,
                           which.clust = NULL,
                           resolution.vals = c(.05, .1, .2, .35),
                           k.val = NULL,
                           do.plot = FALSE,
                           random.seed = 629) {
  # check inputs
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object as input!") }
  # run function
  if (!is.null(which.clust)) {
    reclust_list <- list()
    unique_clusts <- sort(as.integer(unique(seurat.object$seurat_clusters)) - 1)
    # recluster for a single cluster
    if (length(which.clust == 1)) {
      print(sprintf("Identifying subpopulations in cluster %s using %s highly variable genes",
                    which.clust,
                    n.variable.genes))
      temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[1]])
      temp_obj <- SCTransform(temp_obj,
                              vars.to.regress = "percent_MT",
                              seed.use = random.seed,
                              variable.features.n = n.variable.genes,
                              verbose = FALSE)
      temp_obj <- RunPCA(temp_obj,
                         npcs = n.PC,
                         seed.use = 629,
                         verbose = FALSE,
                         features = VariableFeatures(temp_obj))
      if (redo.tSNE) {
        temp_obj <- RunTSNE(temp_obj,
                          reduction = "pca",
                          dims = 1:n.PC,
                          dim.embed = 2,
                          seed.use = random.seed)
      }
      # set k parameter
      if (is.null(k.val)) { k.val <- round(sqrt(ncol(temp_obj))) }
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                k.param = k.val)
      # iterate over resolution parameters and compute silhouette scores to find best re-clustering
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
      # extract best parameters and save results
      names(sil_scores) <- as.character(resolution.vals)
      if (max(sil_scores) > .25) {
        correct_res <- as.numeric(names(sil_scores[sil_scores == max(sil_scores)]))
        print(sprintf("Clustering cells using resolution = %s, which achieved silhouette score: %s",
                      correct_res,
                      round(max(sil_scores), 4)))
        temp_obj <- FindClusters(temp_obj,
                                 resolution = correct_res,
                                 algorithm = 1,
                                 random.seed = random.seed)
        if (do.plot == TRUE) {
          print(DimPlot(temp_obj, reduction = "tsne") + labs(title = sprintf("Reclustering with resolution = %s", correct_res)))
        }
      } else {
        # replace new object w/ original one, as no subpopulations were found
        print(sprintf("Did not find suffcient evidence of subclusters, as the max silhouette score was: %s", round(max(sil_scores), 4)))
        temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
      }
      reclust_list[[1]] <- temp_obj
      # recluster for multiple clusters
    } else if (length(which.clust) > 1) {
      for (clust in seq(which.clust)) {
        print(sprintf("Identifying subpopulations in cluster %s using %s highly variable genes",
                      which.clust[[clust]],
                      n.variable.genes))
        temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[clust]])
        temp_obj <- SCTransform(temp_obj,
                                vars.to.regress = "percent_MT",
                                seed.use = random.seed,
                                variable.features.n = n.variable.genes,
                                verbose = FALSE)
        temp_obj <- RunPCA(temp_obj,
                           npcs = n.PC,
                           seed.use = 629,
                           verbose = FALSE,
                           features = VariableFeatures(temp_obj))
        if (redo.tSNE) {
          temp_obj <- RunTSNE(temp_obj,
                            reduction = "pca",
                            dims = 1:n.PC,
                            dim.embed = 2,
                            seed.use = random.seed)
        }
        # set k parameter
        if (is.null(k.val)) { k.val <- round(sqrt(ncol(temp_obj))) }
        temp_obj <- FindNeighbors(temp_obj,
                                  reduction = "pca",
                                  k.param = k.val)
        # iterate over resolution parameters and compute silhouette scores to find best re-clustering
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
        # extract best parameters and save results
        names(sil_scores) <- as.character(resolution.vals)
        if (max(sil_scores) > .25) {
          correct_res <- as.numeric(names(sil_scores[sil_scores == max(sil_scores)]))
          print(sprintf("Clustering cells using resolution = %s, which achieved silhouette score: %s",
                        correct_res,
                        round(max(sil_scores), 4)))
          temp_obj <- FindClusters(temp_obj,
                                   resolution = correct_res,
                                   algorithm = 1,
                                   random.seed = random.seed)
          if (do.plot == TRUE) {
            print(DimPlot(temp_obj, reduction = "tsne") + labs(title = sprintf("Reclustering with resolution = %s", correct_res)))
          }
        } else {
          # replace new object w/ original one, as no subpopulations were found
          print(sprintf("Did not find suffcient evidence of subclusters, as the max silhouette score was: %s", round(max(sil_scores), 4)))
          temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[clust]])
        }
        reclust_list[[clust]] <- temp_obj
      }
    }
    # prepare results
    names(reclust_list) <- as.character(unlist(which.clust))
  } else { stop("Please provide a list of clusters to analyze.") }

  return(reclust_list)
}
