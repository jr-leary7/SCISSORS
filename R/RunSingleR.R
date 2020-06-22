#' A function to run `SingleR` on a `Seurat` object.
#'
#' This function runs the `SingleR` cell type identification algorithm, with a `Seurat` object as input.
#' @import Seurat
#' @import SingleR
#' @param seurat.object The `Seurat` object containing the cells you'd like to ID.
#' @param ref.data (Optional) The reference dataset whose labels you'd like to assign to your cells. Defaults to the Human Primary Cell Atlas.
#' @param dif.exp.method The differential expression testing method used by `SingleR`. Defaults to the Wilcoxon ranked-sum test.
#' @export
#' @examples
#' callSingleR(seurat.object = pbmc)
#' callSingleR(seurat.object = pbmc, ref.data = my_ref, dif.exp.method = "t")
#' @references
#' Aran *et al* (2019). Reference-based analysis of lung single-cell sequencing data reveals a transitional profibrotic macrophage.

RunSingleR <- function(seurat.object = NULL,
                        ref.data = NULL,
                        dif.exp.method = "wilcox") {
  # test input
  if (is.null(seurat.object@assays$SCT)) {
    stop("Please normalize your data using SCTransform before running SingleR")
  }

  # setup inputs
  norm_counts <- data.frame(seurat.object@assays$SCT@data)
  if (is.null(ref.data)) {
    bulk_ref <- HumanPrimaryCellAtlasData()
    results <- SingleR(test = norm_counts,
                       ref = bulk_ref,
                       labels = bulk_ref$label.main,
                       clusters = seurat.object$seurat_clusters,
                       method = "cluster",
                       de.method = dif.exp.method)
  table(results$labels)


  }


}
