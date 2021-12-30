#' Identify subpopulations in single cell clusters.
#'
#' @name ReclusterCells
#' @author Jack Leary
#' @description This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using \code{\link[Seurat]{SCTransform}}. The function returns a list of \code{Seurat} objects, one for each cluster the user wants to investigate.
#' @importFrom Seurat DefaultAssay SCTransform FindNeighbors FindClusters
#' @param seurat.object The \code{Seurat} object containing cells and their assigned cluster IDs.
#' @param which.clust Which clusters should undergo subpopulation detection analysis? A user-provided list or single integer. Defaults to NULL.
#' @param auto Should the clusters on which to run SCISSORS be determined automatically? If so, \code{which.clust} will be chosen through silhouette score analysis. Not recommended for large datasets as the distance matrix calculation is computationally expensive. Defaults to FALSE.
#' @param merge.clusters (Optional) If multiple clusters are specified, should the clusters be grouped as one before running SCISSORS? Defaults to FALSE.
#' @param n.HVG How many variable genes should be detected in each subcluster? Defaults to 4000.
#' @param n.PC How many PCs should be used as input to non-linear to non-linear dimension reduction and clustering algorithms. Can be provided by the user, or set automatically by \code{\link{ChoosePCs}}. Defaults to "auto".
#' @param redo.embedding (Optional) Should a cluster-specific dimension reduction embeddings be generated? Sometimes subpopulations appear mixed together on the original coordinates, but separate clearly when re-embedded. Defaults to TRUE.
#' @param resolution.vals A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.1, .2, .3, .4).
#' @param k.vals The values of the number of nearest neighbors \emph{k} to be tested. Defaults to c(10, 25, 50).
#' @param cutoff.score The lowest mean silhouette score accepted as evidence of subclusters. Defaults to .25, reasonable values are \[.1, .3\].
#' @param nn.metric (Optional) The distance metric to be used in computing the SNN graph. Defaults to "cosine".
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @seealso \code{\link{ComputeSilhouetteScores}}
#' @export
#' @examples
#' \dontrun{ReclusterCells(seurat.object, which.clust = 5, resolution.vals = c(.1, .2, .5), k.vals = c(10, 20, 30))}
#' \dontrun{ReclusterCells(seurat.object, which.clust = list(0, 3, 5), merge.clusters = TRUE)}

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
                           nn.metric = "cosine",
                           random.seed = 629) {
  # check inputs
  if (is.null(seurat.object) | is.null(which.clust)) { stop("Please provide a Seurat object and clusters to investigate to ReclusterCells().") }
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
    old_which_clust <- which.clust
    which.clust <- 1
  }
  regress_vars <- c()
  if ("percent_MT" %in% colnames(seurat.object@meta.data)) {
    regress_vars <- c(regress_vars, "percent_MT")
  }
  if ("S.Score" %in% colnames(seurat.object@meta.data) && "G2M.Score" %in% colnames(seurat.object@meta.data)) {
    regress_vars <- c(regress_vars, "S.Score", "G2M.Score")
  }
  # determine which dimension reduction algs to run
  dim_red_algs <- NULL
  if ("tsne" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "tsne")
  } else if ("umap" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "umap")
  } else if ("phate" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "phate")
  }
  # iterate and recluster cells
  for (i in seq_along(which.clust)) {
    if (!merge.clusters) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[i]])
    }
    # reprocess data
    if (Seurat::DefaultAssay(temp_obj) != "integrated") {
      if (length(regress_vars) > 0) {
        temp_obj <- Seurat::SCTransform(temp_obj,
                                        vars.to.regress = regress_vars,
                                        variable.features.n = n.HVG,
                                        seed.use = random.seed,
                                        verbose = FALSE)
      } else {
        temp_obj <- Seurat::SCTransform(temp_obj,
                                        variable.features.n = n.HVG,
                                        seed.use = random.seed,
                                        verbose = FALSE)
      }
    }
    temp_obj <- ReduceDimensions(temp_obj,
                                 n.PC = n.PC,
                                 which.algos = dim_red_algs,
                                 random.seed = random.seed)

    # silhouette score various clusterings to find best results
    sil_scores <- c()
    j <- 1
    for (k in seq_along(k.vals)) {
      for (r in seq_along(resolution.vals)) {
        temp_obj <- Seurat::FindNeighbors(temp_obj,
                                          reduction = "pca",
                                          dims = 1:n.PC,
                                          k.param = k.vals[k],
                                          annoy.metric = nn.metric,
                                          nn.method = "annoy",
                                          verbose = FALSE)
        temp_obj <- Seurat::FindClusters(temp_obj,
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
                      paste(old_which_clust, collapse = ", "),
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
      temp_obj <- Seurat::FindNeighbors(temp_obj,
                                        reduction = "pca",
                                        dims = 1:n.PC,
                                        annoy.metric = nn.metric,
                                        nn.method = "annoy",
                                        k.param = best_k,
                                        verbose = FALSE)
      temp_obj <- Seurat::FindClusters(temp_obj,
                                       resolution = best_res,
                                       algorithm = 1,
                                       random.seed = random.seed,
                                       verbose = FALSE)
    } else {
      # replace new object w/ original one, as no subclusters were found
      if (merge.clusters) {
        print(sprintf("Didn't find subclusters in merged clusters %s; max S = %s",
                      paste(old_which_clust, collapse = ", "),
                      round(max(sil_scores), 3)))
        temp_obj <- subset(seurat.object, subset = seurat_clusters %in% which.clust)
      } else {
        print(sprintf("Didn't find subclusters in cluster %s; max S = %s",
                      which.clust[[i]],
                      round(max(sil_scores), 3)))
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
