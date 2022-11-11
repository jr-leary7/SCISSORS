#' Identify subpopulations in single cell clusters.
#'
#' @name ReclusterCells
#' @author Jack Leary
#' @description This function identifies subclusters of cell types by recalculating the *n* most highly variable genes for each cluster using \code{\link[Seurat]{SCTransform}}. The function returns a list of \code{Seurat} objects, one for each cluster the user wants to investigate.
#' @importFrom Seurat DefaultAssay SplitObject SCTransform FindVariableFeatures NormalizeData SelectIntegrationFeatures FindIntegrationAnchors PrepSCTIntegration IntegrateData ScaleData FindNeighbors FindClusters
#' @importFrom future plan
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterEvalQ
#' @param seurat.object The \code{Seurat} object containing cells and their assigned cluster IDs.
#' @param which.clust Which clusters should undergo subpopulation detection analysis? A user-provided list or single integer. Leave NULL if setting \code{auto = TRUE}. Defaults to NULL.
#' @param auto Should the clusters on which to run SCISSORS be determined automatically? If so, \code{which.clust} will be chosen through silhouette score analysis. Not recommended for large datasets as the distance matrix calculation is computationally expensive. Defaults to FALSE.
#' @param merge.clusters If multiple clusters are specified, should the clusters be grouped as one before running SCISSORS? Defaults to FALSE.
#' @param use.parallel Should the \code{Seurat} data reprocessing & the main reclustering loop be parallelized? Defaults to TRUE.
#' @param n.cores The number of cores to be used in parallel computation is \code{use.parallel = TRUE}. Defaults to 3.
#' @param use.sct Should \code{SCTransform} be used for normalization / HVG selection? Defaults to TRUE, otherwise typical log-normalization is used.
#' @param n.HVG How many variable genes should be detected in each subcluster? Defaults to 4000.
#' @param n.PC How many PCs should be used as input to non-linear to non-linear dimension reduction and clustering algorithms. Can be provided by the user, or set automatically by \code{\link{ChoosePCs}}. Defaults to "auto".
#' @param redo.embedding (Optional) Should a cluster-specific dimension reduction embeddings be generated? Sometimes subpopulations appear mixed together on the original coordinates, but separate clearly when re-embedded. Defaults to TRUE.
#' @param resolution.vals A user-defined vector of resolution values to compare when clustering cells. Defaults to c(.1, .2, .3, .4).
#' @param k.vals The values of the number of nearest neighbors \emph{k} to be tested. Defaults to c(10, 25, 50).
#' @param is.integrated Do the data come from multiple samples & need to be re-integrated? See https://github.com/satijalab/seurat/issues/1883 for discussion on this topic. Defaults to FALSE.
#' @param integration.ident If the data are to be re-integrated, what metadata column contains the sample identity? Defaults to NULL.
#' @param cutoff.score The lowest mean silhouette score accepted as evidence of subclusters. Defaults to .25, reasonable values are \[.1, .3\].
#' @param nn.metric (Optional) The distance metric to be used in computing the SNN graph. Defaults to "cosine".
#' @param regress.mt (Optional) Should the percentage of mitochondrial DNA be computed and regressed out? Works for mouse / human gene names. Defaults to FALSE.
#' @param regress.cc (Optional) Should cell cycle scores be computed & regressed out? NOTE: uses human cell cycle genes. Defaults to FALSE.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 312
#' @seealso \code{\link{ComputeSilhouetteScores}}
#' @export
#' @examples
#' \dontrun{ReclusterCells(seurat.object, which.clust = 5, resolution.vals = c(.1, .2, .5), k.vals = c(10, 20, 30))}
#' \dontrun{ReclusterCells(seurat.object, which.clust = list(0, 3, 5), merge.clusters = TRUE)}
#' \dontrun{ReclusterCells(seurat.object, which.clust = list(0, 2), use.parallel = TRUE, n.cores = 6)}

ReclusterCells <- function(seurat.object = NULL,
                           which.clust = NULL,
                           auto = FALSE,
                           merge.clusters = FALSE,
                           use.parallel = TRUE,
                           n.cores = 3,
                           use.sct = TRUE,
                           n.HVG = 4000,
                           n.PC = "auto",
                           redo.embedding = TRUE,
                           resolution.vals = c(.1, .2, .3, .4),
                           k.vals = c(10, 25, 50),
                           is.integrated = FALSE,
                           integration.ident = NULL,
                           cutoff.score = .25,
                           nn.metric = "cosine",
                           regress.mt = FALSE,
                           regress.cc = FALSE,
                           random.seed = 312) {
  # check inputs
  if (is.null(seurat.object)) { stop("Please provide a Seurat object to ReclusterCells().") }
  if (is.null(which.clust) && !auto) { stop("Please provide a vector of clusters or set auto = TRUE.") }
  if (is.integrated && !integration.ident %in% colnames(seurat.object@meta.data)) { stop("integration.ident must exist in Seurat object metadata.") }

  # auto-choose clusters based on silhouette scores to investigate if desired
  if (auto) {
    print("Choosing reclustering candidates automatically.")
    which.clust <- as.integer(names(which(ComputeSilhouetteScores(seurat.object) < cutoff.score)))
  }

  # set up variables to regress out
  regression_vars <- c()
  if (regress.cc) {
    regression_vars <- c(regression_vars, "CC_difference")
  }
  if (regress.mt) {
    regression_vars <- c(regression_vars, "percent_MT")
  }

  # determine which dimension reduction algos to run
  dim_red_algs <- NULL
  if ("tsne" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "tsne")
  } else if ("umap" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "umap")
  } else if ("phate" %in% names(seurat.object@reductions)) {
    dim_red_algs <- c(dim_red_algs, "phate")
  }

  # set up result list, account for case when clusters are to be merged
  reclust_list <- list()
  if (merge.clusters) {
    temp_obj <- subset(seurat.object, subset = seurat_clusters %in% which.clust)
    old_which_clust <- which.clust
    which.clust <- 1L
  }
  # iterate and recluster cells
  for (i in seq_along(which.clust)) {
    if (!merge.clusters) {
      temp_obj <- subset(seurat.object, subset = seurat_clusters == which.clust[[i]])
    }

    # parallelize Seurat pre-processing with future
    if (use.parallel) {
      future::plan("multisession", workers = n.cores)
    }

    # reprocess data, accounting for integration if necessary (otherwise batch effects are usually pretty large & will define clustering)
    if (is.integrated) {
      obj_list <- Seurat::SplitObject(seurat.object, split.by = integration.ident)
      if (use.sct) {
        # SCTransform integration
        obj_list <- purrr::map(obj_list, SCTransform)
        int_features <- Seurat::SelectIntegrationFeatures(obj_list,
                                                          nfeatures = n.HVG,
                                                          verbose = FALSE)
        obj_list <- Seurat::PrepSCTIntegration(obj_list,
                                               anchor.features = int_features,
                                               verbose = FALSE)
        int_anchors <- Seurat::FindIntegrationAnchors(obj_list,
                                                      normalization.method = "SCT",
                                                      anchor.features = int_features,
                                                      verbose = FALSE)
        temp_obj <- Seurat::IntegrateData(anchorset = int_anchors,
                                          normalization.method = "SCT",
                                          verbose = FALSE)
        # not sure really how to regress out effects w/ SCT after integration, so using this
        if (length(regression_vars) > 0) {
          temp_obj <- Seurat::ScaleData(temp_obj,
                                        vars.to.regress = regression_vars,
                                        model.use = "negbinom",
                                        verbose = FALSE)
        } else {
          temp_obj <- Seurat::ScaleData(temp_obj, verbose = FALSE)
        }
      } else {
        # log-normalization integration
        obj_list <- purrr::map(obj_list, function(x) {
          x %>%
            Seurat::NormalizeData(verbose = FALSE) %>%
            Seurat::FindVariableFeatures(nfeatures = n.HVG, verbose = FALSE)
        })
        int_features <- Seurat::SelectIntegrationFeatures(obj_list,
                                                          nfeatures = n.HVG,
                                                          verbose = FALSE)
        int_anchors <- Seurat::FindIntegrationAnchors(obj_list,
                                                      anchor.features = int_features,
                                                      normalization.method = "LogNormalize",
                                                      verbose = FALSE)
        temp_obj <- Seurat::IntegrateData(anchorset = int_anchors,
                                          normalization.method = "LogNormalize",
                                          verbose = FALSE)
        if (length(regression_vars) > 0) {
          temp_obj <- Seurat::ScaleData(temp_obj,
                                        vars.to.regress = regression_vars,
                                        model.use = "negbinom",
                                        verbose = FALSE)
        } else {
          temp_obj <- Seurat::ScaleData(temp_obj, verbose = FALSE)
        }
      }
    } else {
      if (use.sct) {
        # SCTransform normalization (no integration)
        if (length(regression_vars) > 0) {
          temp_obj <- Seurat::SCTransform(temp_obj,
                                          vars.to.regress = regression_vars,
                                          variable.features.n = n.HVG,
                                          seed.use = random.seed,
                                          verbose = FALSE)
        } else {
          temp_obj <- Seurat::SCTransform(temp_obj,
                                          variable.features.n = n.HVG,
                                          seed.use = random.seed,
                                          verbose = FALSE)
        }
      } else {
        # log-normalization (no integration)
        temp_obj <- Seurat::NormalizeData(temp_obj,
                                          normalization.method = "LogNormalize",
                                          scale.factor = 10000,
                                          verbose = FALSE)
        temp_obj <- Seurat::FindVariableFeatures(temp_obj,
                                                 selection.method = "vst",
                                                 nfeatures = n.HVG,
                                                 verbose = FALSE)
        if (length(regression_vars) > 0) {
          temp_obj <- Seurat::ScaleData(temp_obj,
                                        vars.to.regress = regression_vars,
                                        model.use = "negbinom",
                                        verbose = FALSE)
        } else {
          temp_obj <- Seurat::ScaleData(temp_obj, verbose = FALSE)
        }
      }
    }
    if (redo.embedding) {
      temp_obj <- ReduceDimensions(temp_obj,
                                   n.PC = n.PC,
                                   which.algos = dim_red_algs,
                                   random.seed = random.seed)
    }
    if (use.parallel) {
      future:::ClusterRegistry("stop")
    }
    # silhouette score various clusterings to find best results
    if (use.parallel) {
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
    } else {
      cl <- foreach::registerDoSEQ()
      set.seed(random.seed)
    }

    param_grid <- expand.grid(k.vals, resolution.vals)
    colnames(param_grid) <- c("K", "R")
    # loop over parameter values in parallel
    sil_scores <- foreach::foreach(i = seq(nrow(param_grid)),
                                   .combine = "c",
                                   .multicombine = ifelse(nrow(param_grid) > 1, TRUE, FALSE),
                                   .maxcombine = ifelse(nrow(param_grid) > 1, nrow(param_grid), 2),
                                   .packages = c("SCISSORS", "Seurat"),
                                   .verbose = FALSE) %dopar% {
      k <- param_grid$K[i]
      r <- param_grid$K[i]
      temp_obj <- Seurat::FindNeighbors(temp_obj,
                                        reduction = "pca",
                                        dims = 1:n.PC,
                                        k.param = k,
                                        annoy.metric = nn.metric,
                                        nn.method = "annoy",
                                        verbose = FALSE) %>%
                  Seurat::FindClusters(resolution = r,
                                       random.seed = random.seed,
                                       algorithm = 1,
                                       verbose = FALSE)
      if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
        sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
        mean_sil <- mean(sil_res)
      } else {
        # neutral placeholder value for the case when the number of identified clusters is 1
        mean_sil <- 0
      }
      mean_sil
    }
    names(sil_scores) <- as.character(paste0("k", param_grid$K, "r", param_grid$R))
    # end parallelization & clean up
    sink(tempfile())
    if (use.parallel) {
      parallel::clusterEvalQ(cl, expr = {
        rm(list = ls(all.names = TRUE)); gc(verbose = FALSE, full = TRUE)
      })
      parallel::stopCluster(cl)
    }
    rm(cl)
    gc(verbose = FALSE, full = TRUE)
    sink()

    # OLD LOOP - keeping code here until release is finalized
    # sil_scores <- c()
    # j <- 1
    # for (k in seq_along(k.vals)) {
    #   for (r in seq_along(resolution.vals)) {
    #     temp_obj <- Seurat::FindNeighbors(temp_obj,
    #                                       reduction = "pca",
    #                                       dims = 1:n.PC,
    #                                       k.param = k.vals[k],
    #                                       annoy.metric = nn.metric,
    #                                       nn.method = "annoy",
    #                                       verbose = FALSE)
    #     temp_obj <- Seurat::FindClusters(temp_obj,
    #                                      resolution = resolution.vals[r],
    #                                      random.seed = random.seed,
    #                                      algorithm = 1,
    #                                      verbose = FALSE)
    #     if (length(unique(levels(temp_obj$seurat_clusters))) > 1) {
    #       sil_res <- ComputeSilhouetteScores(seurat.obj = temp_obj)
    #       mean_sil <- mean(sil_res)
    #       sil_scores[j] <- mean_sil
    #     } else {
    #       # neutral placeholder value for the case when the number of identified clusters is 1
    #       sil_scores[j] <- 0
    #     }
    #     names(sil_scores)[j] <- as.character(paste0("k", k.vals[k], "r", resolution.vals[r]))
    #     j <- j + 1
    #   }
    # }

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
