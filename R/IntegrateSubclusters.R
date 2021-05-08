#' A function to integrate subcluster identities into the original Seurat object.
#'
#' This function takes a list of outputs from ReclusterCells() and integrates the new subcluster identities into the original Seurat object.
#' @import Seurat
#' @import dplyr
#' @param original.object The original Seurat object. Defaults to NULL.
#' @param reclust.results A list of reclustering results as output from ReclusterCells(). Defaults to NULL.
#' @param do.plot Should the results be plotted on a dimension reduction plot? Defaults to FALSE.
#' @export
#' @examples
#' IntegrateSubclusters(original.object = pbmc, reclust.results = my_subclusts)

IntegrateSubclusters <- function(original.object = NULL,
                                 reclust.results = NULL,
                                 do.plot = FALSE) {
  # check inputs
  if (any(sapply(c(original.object, reclust.results), is.null))) stop("Please provide a Seurat object and list of reclustering results.")
  if (class(reclust.results) != "list") stop("reclust.results must be a list.")
  if (any(sapply(reclust.results, class) != "Seurat")) stop("All elements of reclust.results must be Seurat objects.")
  max_clust <- max(as.numeric(original.object$seurat_clusters) - 1)
  # identify new subclusters
  cell_df <- NULL
  for (i in seq_along(reclust.results)) {
    n_clust <- length(unique(reclust.results[[i]]$seurat_clusters))
    if (n_clust == 1) {
      next
    } else {
      meta_df <- reclust.results[[i]]@meta.data
      unique_subclust <- sort(unique(as.numeric(meta_df$seurat_clusters)) - 1)
      for (j in seq_along(unique_subclust)) {
        subclust_df <- data.frame(Cell = rownames(meta_df[meta_df$seurat_clusters == unique_subclust[j], ]),
                                  ClustID = max_clust + 1)
        cell_df <- rbind(cell_df, subclust_df)
        max_clust <- max_clust + 1
      }
    }
  }
  # integrate new subclusters
  new_clusts <- unique(cell_df$ClustID)
  original.object@meta.data$seurat_clusters_original <- original.object@meta.data$seurat_clusters
  original.object@meta.data$seurat_clusters <- as.numeric(original.object@meta.data$seurat_clusters_original) - 1
  for (k in seq_along(new_clusts)) {
    clust_cells <- cell_df[cell_df$ClustID == new_clusts[k], ]$Cell
    original.object@meta.data$seurat_clusters <- case_when(rownames(original.object@meta.data) %in% clust_cells ~ new_clusts[k],
                                                           TRUE ~ original.object@meta.data$seurat_clusters)
  }
  clusts_to_fix <- sort(unique(original.object$seurat_clusters))
  for (l in seq_along(clusts_to_fix)) {
    original.object@meta.data[original.object@meta.data$seurat_clusters == clusts_to_fix[l], ]$seurat_clusters <- l
  }
  original.object@meta.data$seurat_clusters <- as.factor(original.object@meta.data$seurat_clusters - 1)
  Idents(original.object) <- "seurat_clusters"

  # plot results if desired
  if (do.plot) {
    p <- DimPlot(original.object) +
         labs(title = "Integrated Subclusters") +
         theme_yehlab() +
         theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }

  return(original.object)
}
