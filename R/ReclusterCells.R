#' Identify subpopulations in single cell clusters.
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using `sctransform` as implemented in `Seurat`. The function returns a list of `Seurat` objects, one for each cluster the user wants to investigate.
#' @import Seurat
#' @importFrom ggplot2 labs
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param which.clust Which clusters should undergo subpopulation detection analysis? A user-provided list.
#' @param merge.clusters (Optional). If multiple clusters are specified, should the clusters be re-clustered as one? Defaults to FALSE.
#' @param n.variable.genes How many variable genes should be detected in each subcluster? Defaults to 4000.
#' @param n.PC How many PCs should be used as input to non-linear to non-linear dimension reduction and clustering algorithms. Defaults to 10.
#' @param redo.embedding (Optional) Should a cluster-specific dimension reduction embeddings be generated? Sometimes subpopulations appear mixed together on the original coordinates, but separate clearly when re-embedded. Defaults to TRUE.
#' @param which.dim.reduc (Optional). Which non-linear dimension reduction algorithms should be used? Supports "tsne", "umap", "phate", and "all". Plots will be generated using the t-SNE embedding. Defaults to c("tsne", "umap"), as most users will likely not have `phateR` installed.
#' @param resolution.vals (Optional) A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.05, .1, .15, .2, .35).
#' @param k.vals (Optional) The parameters *k* to be tested .
#' @param do.plot (Optional) Should t-SNE plots of the various reclusterings be plotted for visual inspection by the user? Defaults to FALSE.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' ReclusterCells(seurat.object, resolution.vals = c(.1, .2, .5))
#' ReclusterCells(seurat.object, which.dim.reduc = c("tsne", "phate"))
#' ReclusterCells(seurat.object, which.clust = list(0, 3, 5), do.plot = TRUE)

ReclusterCells <- function(seurat.object = NULL,
                           which.clust = NULL,
                           merge.clusters = FALSE,
                           n.variable.genes = 4000,
                           n.PC = 10,
                           which.dim.reduc = c("tsne", "umap"),
                           redo.embedding = TRUE,
                           resolution.vals = c(.05, .1, .2, .35),
                           k.vals = c(10, 25, 50),
                           do.plot = FALSE,
                           random.seed = 629) {
  # check inputs
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object as input!") }
  # run function
  if (!is.null(which.clust)) {
    reclust_list <- list()
    # recluster for a single cluster or multiple clusters merged into one
    if (length(which.clust) == 1 | merge.clusters == TRUE) {
      # create subset
      if (length(which.clust) == 1) {
        temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[1]])
      } else if (merge.clusters == TRUE) {
        temp_obj <- subset(seurat.object, subset = seurat_clusters %in% which.clust)
      }
      temp_obj <- SCTransform(temp_obj,
                              vars.to.regress = "percent_MT",
                              seed.use = random.seed,
                              variable.features.n = n.variable.genes,
                              verbose = FALSE)
      temp_obj <- RunPCA(temp_obj,
                         npcs = n.PC,
                         features = VariableFeatures(temp_obj),
                         seed.use = random.seed,
                         verbose = FALSE)
      # re-run nonlinear dimension reduction
      if (redo.embedding) {
        if (which.dim.reduc == "all") {
          temp_obj <- RunTSNE(temp_obj,
                              reduction = "pca",
                              dims = 1:n.PC,
                              dim.embed = 2,
                              seed.use = random.seed)
          temp_obj <- RunUMAP(temp_obj,
                              reduction = "pca",
                              dims = 1:n.PC,
                              umap.method = "uwot",
                              n.components = 2,
                              metric = "cosine",
                              seed.use = random.seed,
                              verbose = FALSE)
          pca_df <- data.frame(Embeddings(temp_obj, reduction = "pca"))
          phate_res <- phate(pca_df,
                             ndim = 2,
                             mds.solver = "smacof",
                             knn.dist.method = "cosine",
                             mds.dist.method = "cosine",
                             npca = NULL,
                             seed = random.seed,
                             verbose = FALSE)
          phate_obj <- CreateDimReducObject(embeddings = phate_res$embedding,
                                            assay = "SCT",
                                            key = "PHATE_",
                                            global = TRUE)
          temp_obj@reductions$phate <- phate_obj
        }
        if ("tsne" %in% which.dim.reduc) {
          temp_obj <- RunTSNE(temp_obj,
                              reduction = "pca",
                              dims = 1:n.PC,
                              dim.embed = 2,
                              seed.use = random.seed)
        }
        if ("umap" %in% which.dim.reduc) {
          temp_obj <- RunUMAP(temp_obj,
                              reduction = "pca",
                              dims = 1:n.PC,
                              umap.method = "uwot",
                              n.components = 2,
                              metric = "cosine",
                              seed.use = random.seed,
                              verbose = FALSE)
        }
        if ("phate" %in% which.dim.reduc) {
          pca_df <- data.frame(Embeddings(temp_obj, reduction = "pca"))
          phate_res <- phate(pca_df,
                             ndim = 2,
                             mds.solver = "smacof",
                             knn.dist.method = "cosine",
                             mds.dist.method = "cosine",
                             npca = NULL,
                             seed = random.seed,
                             verbose = FALSE)
          phate_obj <- CreateDimReducObject(embeddings = phate_res$embedding,
                                            assay = "SCT",
                                            key = "PHATE_",
                                            global = TRUE)
          temp_obj@reductions$phate <- phate_obj
        }
      }
      # set max k parameter
      k_max <- round(sqrt(ncol(temp_obj)))
      if (k_max > max(k.vals)) {
        k.vals <- c(k.vals, k_max)
      }
      # iterate over resolution parameters and compute silhouette scores to find best re-clustering
      sil_scores <- c()
      i <- 1
      for (k in seq(k.vals)) {
        for (r in seq(resolution.vals)) {
          temp_obj <- FindNeighbors(temp_obj,
                                    reduction = "pca",
                                    k.param = k.vals[k],
                                    verbose = FALSE)
          temp_obj <- FindClusters(temp_obj,
                                   resolution = resolution.vals[r],
                                   random.seed = random.seed,
                                   algorithm = 1,
                                   verbose = FALSE)
          if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
            sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
            mean_sil <- mean(sil_res)
            sil_scores[i] <- mean_sil
          } else {
            # neutral placeholder value for the case when the number of identified clusters is 1
            sil_scores[i] <- 0
          }
          names(sil_scores)[i] <- as.character(paste0("k", k.vals[k], "r", resolution.vals[r]))
          i <- i + 1
        }
      }
      # extract best parameters and save results
      if (max(sil_scores) > .25) {
        best_params <- names(sil_scores[sil_scores == max(sil_scores)])
        if (length(best_params) > 1) {
          best_param_index <- ceiling(length(best_params) / 2)
          best_k <- strsplit(best_params[best_param_index], "r")[[1]][1]
          best_k <- as.numeric(strsplit(best_k, "k")[[1]][2])
          best_res <- as.numeric(strsplit(best_params[best_param_index], "r")[[1]][2])
        } else {
          best_k <- strsplit(best_params, "r")[[1]][1]
          best_k <- as.numeric(strsplit(best_k, "k")[[1]][2])
          best_res <- as.numeric(strsplit(best_params, "r")[[1]][2])
        }
        # cluster cells using best parameters
        print(sprintf("Reclustering cells in cluster %s using k = %s & resolution = %s, which achieved silhouette score: %s",
                      which.clust[[1]],
                      best_k,
                      best_res,
                      round(max(sil_scores), 3)))
        temp_obj <- FindNeighbors(temp_obj,
                                  reduction = "pca",
                                  k.param = best_k,
                                  verbose = FALSE)
        temp_obj <- FindClusters(temp_obj,
                                 resolution = best_res,
                                 algorithm = 1,
                                 random.seed = random.seed,
                                 verbose = FALSE)
        if (do.plot) {
          print(DimPlot(temp_obj, reduction = "tsne") + labs(title = sprintf("Reclustering of cluster %s using k = %s & resolution = %s",
                                                                             which.clust[[1]],
                                                                             best_k,
                                                                             best_res)))
        }
      } else {
        # replace new object w/ original one, as no subpopulations were found
        print(sprintf("Did not find suffcient evidence of subclusters in cluster %s, as the max silhouette score was: %s",
                      which.clust[[1]],
                      round(max(sil_scores), 3)))
        temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[1]])
      }
      reclust_list[[1]] <- temp_obj

      # recluster for multiple clusters
    } else if (length(which.clust) > 1 & merge.clusters == FALSE) {
      for (clust in seq(which.clust)) {
        temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[clust]])
        temp_obj <- SCTransform(temp_obj,
                                vars.to.regress = "percent_MT",
                                seed.use = random.seed,
                                variable.features.n = n.variable.genes,
                                verbose = FALSE)
        temp_obj <- RunPCA(temp_obj,
                           npcs = n.PC,
                           features = VariableFeatures(temp_obj),
                           seed.use = 629,
                           verbose = FALSE)
        # re-run nonlinear dimension reduction
        if (redo.embedding) {
          if (which.dim.reduc == "all") {
            temp_obj <- RunTSNE(temp_obj,
                                reduction = "pca",
                                dims = 1:n.PC,
                                dim.embed = 2,
                                seed.use = random.seed)
            temp_obj <- RunUMAP(temp_obj,
                                reduction = "pca",
                                umap.method = "uwot",
                                n.components = 2,
                                metric = "cosine",
                                seed.use = random.seed,
                                verbose = FALSE)
            pca_df <- data.frame(Embeddings(temp_obj, reduction = "pca"))
            phate_res <- phate(pca_df,
                               ndim = 2,
                               mds.solver = "smacof",
                               knn.dist.method = "cosine",
                               mds.dist.method = "cosine",
                               npca = NULL,
                               seed = random.seed,
                               verbose = FALSE)
            phate_obj <- CreateDimReducObject(embeddings = phate_res$embedding,
                                              assay = "SCT",
                                              key = "PHATE_",
                                              global = TRUE)
            temp_obj@reductions$phate <- phate_obj
          }
          if ("tsne" %in% which.dim.reduc) {
            temp_obj <- RunTSNE(temp_obj,
                                reduction = "pca",
                                dims = 1:n.PC,
                                dim.embed = 2,
                                seed.use = random.seed)
          }
          if ("umap" %in% which.dim.reduc) {
            temp_obj <- RunUMAP(temp_obj,
                                reduction = "pca",
                                dims = 1:n.PC,
                                umap.method = "uwot",
                                n.components = 2,
                                metric = "cosine",
                                seed.use = random.seed,
                                verbose = FALSE)
          }
          if ("phate" %in% which.dim.reduc) {
            pca_df <- data.frame(Embeddings(temp_obj, reduction = "pca"))
            phate_res <- phate(pca_df,
                               ndim = 2,
                               mds.solver = "smacof",
                               knn.dist.method = "cosine",
                               mds.dist.method = "cosine",
                               npca = NULL,
                               seed = random.seed,
                               verbose = FALSE)
            phate_obj <- CreateDimReducObject(embeddings = phate_res$embedding,
                                              assay = "SCT",
                                              key = "PHATE_",
                                              global = TRUE)
            temp_obj@reductions$phate <- phate_obj
          }
        }
        # set max k parameter
        k_max <- round(sqrt(ncol(temp_obj)))
        if (k_max > max(k.vals)) {
          k.vals <- c(k.vals, k_max)
        }
        # iterate over resolution parameters and compute silhouette scores to find best re-clustering
        sil_scores <- c()
        i <- 1
        for (k in seq(k.vals)) {
          for (r in seq(resolution.vals)) {
            temp_obj <- FindNeighbors(temp_obj,
                                      reduction = "pca",
                                      k.param = k.vals[k],
                                      verbose = FALSE)
            temp_obj <- FindClusters(temp_obj,
                                     resolution = resolution.vals[r],
                                     random.seed = random.seed,
                                     algorithm = 1,
                                     verbose = FALSE)
            if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
              sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
              mean_sil <- mean(sil_res)
              sil_scores[i] <- mean_sil
            } else {
              # neutral placeholder value for the case when the number of identified clusters is 1
              sil_scores[i] <- 0
            }
            names(sil_scores)[i] <- as.character(paste0("k", k.vals[k], "r", resolution.vals[r]))
            i <- i + 1
          }
        }
        # extract best parameters and save results
        if (max(sil_scores) > .25) {
          best_params <- names(sil_scores[sil_scores == max(sil_scores)])
          if (length(best_params) > 1) {
            best_param_index <- ceiling(length(best_params) / 2)
            best_k <- strsplit(best_params[best_param_index], "r")[[1]][1]
            best_k <- as.numeric(strsplit(best_k, "k")[[1]][2])
            best_res <- as.numeric(strsplit(best_params[best_param_index], "r")[[1]][2])
          } else {
            best_k <- strsplit(best_params, "r")[[1]][1]
            best_k <- as.numeric(strsplit(best_k, "k")[[1]][2])
            best_res <- as.numeric(strsplit(best_params, "r")[[1]][2])
          }
          # cluster cells using best parameters
          print(sprintf("Reclustering cells in cluster %s using k = %s & resolution = %s, which achieved silhouette score: %s",
                        which.clust[[clust]],
                        best_k,
                        best_res,
                        round(max(sil_scores), 3)))
          temp_obj <- FindNeighbors(temp_obj,
                                    reduction = "pca",
                                    k.param = best_k,
                                    verbose = FALSE)
          temp_obj <- FindClusters(temp_obj,
                                   resolution = best_res,
                                   algorithm = 1,
                                   random.seed = random.seed,
                                   verbose = FALSE)
          if (do.plot) {
            print(DimPlot(temp_obj, reduction = "tsne") + labs(title = sprintf("Reclustering of cluster %s using k = %s & resolution = %s",
                                                                               which.clust[[clust]],
                                                                               best_k,
                                                                               best_res)))
          }
        } else {
          # replace new object w/ original one, as no subpopulations were found
          print(sprintf("Did not find suffcient evidence of subclusters in cluster %s, as the max silhouette score was: %s",
                        which.clust[[clust]],
                        round(max(sil_scores), 3)))
          temp_obj <- subset(seurat.object, subset = seurat_clusters == clust)
        }
        reclust_list[[clust]] <- temp_obj
      }
    } else { stop("Please provide a list of clusters to analyze.") }
    names(reclust_list) <- as.character(unlist(which.clust))
  }
  return(reclust_list)
}
