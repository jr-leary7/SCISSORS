#' Identify subpopulations in single cell clusters.
#'
#' This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using `sctransform` as implemented in `Seurat`. The function returns a list of `Seurat` objects, one for each cluster the user wants to investigate.
#' @import Seurat
#' @param seurat.object The `Seurat` object containing cells and their assigned cluster IDs.
#' @param which.clust Which clusters should undergo subpopulation detection analysis? A user-provided list or single integer. Defaults to NULL.
#' @param auto Should the clusters on which to run SCISSORS be determined automatically? If so, `which.clust` will be chosen through silhouette score analysis. Not recommended for large datasets as the distance matrix calculation is computationally expensive. Defaults to FALSE.
#' @param merge.clusters (Optional). If multiple clusters are specified, should the clusters be grouped as one before running SCISSORS? Defaults to FALSE.
#' @param n.HVG How many variable genes should be detected in each subcluster? Defaults to 4000.
#' @param n.PC How many PCs should be used as input to non-linear to non-linear dimension reduction and clustering algorithms. Can be provided by the user, or set automatically by `ChoosePCs()`. Defaults to "auto".
#' @param redo.embedding (Optional) Should a cluster-specific dimension reduction embeddings be generated? Sometimes subpopulations appear mixed together on the original coordinates, but separate clearly when re-embedded. Defaults to TRUE.
#' @param resolution.vals (Optional) A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.1, .2, .3, .4).
#' @param k.vals (Optional) The parameters *k* to be tested. Defaults to c(10, 25, 50).
#' @param cutoff.score (Optional) The lowest mean silhouette score accepted as evidence of subclusters. Defaults to .25, reasonable values are [2, .3].
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' ReclusterCells(seurat.object, which.clust = 5, resolution.vals = c(.1, .2, .5), k.vals = c(10, 20, 30))
#' ReclusterCells(seurat.object, which.clust = list(0, 3, 5), merge.clusters = TRUE

ReclusterCells <- function(seurat.object = NULL,
                           which.clust = NULL,
                           auto = FALSE,
                           merge.clusters = FALSE,
                           n.HVG = 4000,
                           n.PC = "auto",
                           redo.embedding = TRUE,
                           resolution.vals = c(.1, .2, .3, .4),
                           k.vals = c(10, 25, 50),
                           cutoff.score = .25,
                           random.seed = 629) {
  # check inputs
  if (any(sapply(c(seurat.object, which.clust), is.null))) stop("Please provide a Seurat object and clusters to investigate.")

  # auto-choose clusters to investigate if desired
  if (auto) {
    print("Choosing clusters automatically.")
    scores <- ComputeSilhouetteScores(seurat.object)
    which.clust <- which(scores < .5)
  }

  # set up result list, account for case when clusters are to be merged, identify covariates
  reclust_list <- list()
  if (merge.clusters) {
    temp_obj <- subset(seurat.object, subset = seurat_clusters %in% which.clust)
    which.clust <- 1
  }
  regress_vars <- c()
  regress_vars <- ifelse("percent_MT" %in% colnames(seurat.object@meta.data),
                         c(regress_vars, "percent_MT"),
                         regress_vars)
  regress_vars <- ifelse("S.Score" %in% colnames(seurat.object@meta.data) && "G2M.Score" %in% colnames(seurat.object@meta.data),
                         c(regress_vars, "S.Score", "G2M.Score"),
                         regress_vars)
  dim_red_algs <- NULL
  dim_red_algs <- ifelse("tsne" %in% names(seurat.object@reductions), c(dim_red_algs, "tsne"), dim_red_algs)
  dim_red_algs <- ifelse("umap" %in% names(seurat.object@reductions), c(dim_red_algs, "umap"), dim_red_algs)
  dim_red_algs <- ifelse("phate" %in% names(seurat.object@reductions), c(dim_red_algs, "phate"), dim_red_algs)
  # iterate and recluster cells
  for (i in seq_along(which.clust)) {
    if (!merge.clusters) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[i]])
    }
    # reprocess data
    temp_obj <- SCTransform(temp_obj,
                            vars.to.regress = regress_vars,
                            seed.use = random.seed,
                            variable.features.n = n.HVG,
                            verbose = FALSE)
    temp_obj <- ReduceDimensions(temp_obj,
                                 n.PC = n.PC,
                                 which.algs = dim_red_algs,
                                 seed = random.seed)

    # silhouette score various clusterings to find best results
    sil_scores <- c()
    j <- 1
    for (k in seq_along(k.vals)) {
      for (r in seq_along(resolution.vals)) {
        temp_obj <- FindNeighbors(temp_obj,
                                  reduction = "pca",
                                  dims = 1:n.PC,
                                  k.param = k.vals[k],
                                  annoy.metric = "cosine",
                                  nn.method = "annoy",
                                  verbose = FALSE)
        temp_obj <- FindClusters(temp_obj,
                                 resolution = resolution.vals[r],
                                 random.seed = random.seed,
                                 algorithm = 1,
                                 verbose = FALSE)
        if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
          sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
          mean_sil <- mean(sil_res)
          sil_scores[j] <- mean_sil
        } else {
          # neutral placeholder value for the case when the number of identified clusters is 1
          sil_scores[j] <- 0
        }
        names(sil_scores)[j] <- as.character(paste0("k", k.vals[k], "r", resolution.vals[r]))
        j <- j + 1
      }
    }

    # extract best parameter set
    if (max(sil_scores) > cutoff.score) {
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
      # cluster cells using best parameter set
      if (merge.clusters) {
        print(sprintf("Reclustering cells in clusters %s using k = %s & resolution = %s; S = %s",
                      paste(which.clust, collapse = ", "),
                      best_k,
                      best_res,
                      round(max(sil_scores), 3)))
      } else {
        print(sprintf("Reclustering cells in cluster %s using k = %s & resolution = %s; S = %s",
                      which.clust[[1]],
                      best_k,
                      best_res,
                      round(max(sil_scores), 3)))
      }
      temp_obj <- FindNeighbors(temp_obj,
                                reduction = "pca",
                                dims = 1:n.PC,
                                annoy.metric = "cosine",
                                nn.method = "annoy",
                                k.param = best_k,
                                verbose = FALSE)
      temp_obj <- FindClusters(temp_obj,
                               resolution = best_res,
                               algorithm = 1,
                               random.seed = random.seed,
                               verbose = FALSE)
    } else {
      # replace new object w/ original one, as no subclusters were found
      if (merge.clusters) {
        print(sprintf("Didn't find subclusters in merged clusters %s; max S = %s"),
              paste(which.clust, collapse = ", "),
              round(max(sil_scores), 3))
        temp_obj <- subset(seurat.object, subset = seurat_clusters %in% which.clust)
      } else {
        print(sprintf("Didn't find subclusters in cluster %s; max S = %s"),
              which.clust[[i]],
              round(max(sil_scores), 3))
        temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[i]])
      }
    }
    reclust_list[[i]] <- temp_obj
  }

  # add names to list of Seurat objects
  if (!merge.clusters && length(which.clust) > 1) {
    names(reclust_list) <- as.character(unlist(which.clust))
  }
  # return single Seurat object or list of objects
  if (length(which.clust) == 1 || merge.clusters) {
    val <- reclust_list[[1]]
  } else {
    val <- reclust_list
  }
  return(val)
}
