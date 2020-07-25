#' Identify marker genes for previosuly identified subpopulations.
#'
#' This function determines which cells characterize the subpopulations identified using `ReclusterCells`. It is intended to be run on a single re-clustered `Seurat` object at a time, though if you wish you could
#' iterate over the list of reclustering results, and save the outputs from this function in a matching array of lists. The function returns a list of dataframes, one dataframe per cluster, containing normal and Bonferroni-adjust
#' p-values, gene prevalence, and effect size in the form of log2 fold change.
#' @import Seurat
#' @param seurat.object The original `Seurat` object containing the entire cell population and related metadata.
#' @param reclust.data The list of reclustered `Seurat` objects as returned by `ReclusterCells`.
#' @param which.clust The cluster (with subpopulations) for which you'd like to obtain marker genes with respect to the entire cell population.
#' @param diff.exp.test The test used to calculate differential expression using `FindMarkers`. Defaults to "wilcox".
#' @param logfc.thresh The lo2 fold-change cutoff used when performing differential expression analysis. Defaults to 2.
#' @param random.seed (Optional) The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' FindSubpopulationMarkers(seurat.object = pbmc_small, reclust.data = reclust_results)
#' FindSubpopulationMarkers(seurat.object = pbmc_small, reclust.data = reclust_results, diff.exp.test = "t")

FindSubpopulationMarkers <- function(seurat.object = NULL,
                                     reclust.data = NULL,
                                     which.clust = NULL,
                                     diff.exp.test = "wilcox",
                                     logfc.thresh = 2,
                                     random.seed = 629) {
  # check inputs
  if (is.null(seurat.object) | is.null(reclust.data) | is.null(which.clust)) { stop("Please provide the correct inputs.") }
  # run function
  temp_obj <- reclust.data[[1]]
  marker_gene_list <- list()
  unique_clusts <- sort(as.integer(unique(temp_obj$seurat_clusters)) - 1)
  # calculate marker genes for each subpopulation cluster vs. all other cells
  for (i in seq(unique_clusts)) {
    clust_cells <- rownames(temp_obj@meta.data[temp_obj@meta.data$seurat_clusters == unique_clusts[i], ])
    big_temp_obj <- seurat.object
    big_temp_obj@meta.data$clust_indicator <- ifelse(rownames(big_temp_obj@meta.data) %in% clust_cells, "Subpopulation", "Other")
    print(sprintf("Finding markers genes for subcluster %s using the %s test",
                  unique_clusts[i],
                  diff.exp.test))
    markers <- FindMarkers(big_temp_obj,
                           slot = "data",
                           assay = "SCT",
                           ident.1 = "Subpopulation",
                           ident.2 = "Other",
                           group.by = "clust_indicator",
                           logfc.threshold = logfc.thresh,
                           verbose = FALSE,
                           test.use = "wilcox",
                           random.seed = random.seed)
    marker_gene_list[[i]] <- markers
  }
  names(marker_gene_list) <- as.character(unique_clusts)

  return(marker_gene_list)
}
