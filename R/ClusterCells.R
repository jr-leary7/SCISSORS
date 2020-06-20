#' Iteratively recluster cells.
#'
#' This function iteratively reclusters cells, each time recalculating the highly variable genes that define each cell population. The goal is to provide more accurate cell subtype recognition, as the highly variable genes that define broad cell types most likely do not correctly represent cell subtypes and rare cell populations.
#' @import Seurat
#' @import SingleR
#' @import SummarizedExperiment
#' @param seurat.object The Seurat object containing the cells you'd like to analyze.
#' @param n.variable.genes The number of variable genes to find at each step. Defaults to 3000.
#' @param initial.resolution The initial resolution parameter used in the `FindClusters` function. Defaults to 0.5.
#' @param perform.SingleR Should broad cell types be determined using the `SingleR` package? Defaults to TRUE.
#' @param use.SingleR.clusters If cell types are determined by `SingleR`, should those types be used as clusters? Defaults to TRUE.
#' @param species If cell types are to be determined using `SingleR`, what species are the cells? Defaults to "human".
#' @param ref.data If cell types are to be determined using `SingleR`, a user-provided dataset can be used as a reference. Defaults to the data from Mabbott *et al* (2013): the Human Primary Cell Atlas.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' ClusterCells(seurat.object = pbmc, n.variable.genes = 5000, initial.resolution = .75, ref.data = cell_ref)
#' ClusterCells(seurat.object = baron_pancreas, initial.resolution = .5, perform.SingleR = TRUE, species = "mouse")
#' @references
#' Stuart *et al* (2019). Comprehensive integration of single-cell data. *Cell*.
#'
#' Aran *et al* (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*.

ClusterCells <- function(seurat.object = NULL,
                         n.variable.genes = 3000,
                         initial.resolution = .5,
                         perform.SingleR = TRUE,
                         species = "human",
                         ref.data = NULL,
                         random.seed = 629) {
  # check arguments & assays present in Seurat object
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object!") }
  if (is.null(seurat.object@assays$SCT) & length(VariableFeatures(seurat.object)) == 0) {
    print("Normalizing data, scaling, and selecting variable features using SCTransform")

    # check if % mito DNA exists in Seurat object & regress out if so
    if (any(grepl("MT|mt", colnames(seurat.object@meta.data)))) {
      col_loc <- which(grepl("MT|mt", colnames(seurat.object@meta.data)))
      col_name <- colnames(seurat.object$meta.data)[col_loc]
      print("Regressing out % mitochondrial DNA")
      seurat.object <- SCTransform(seurat.object,
                                   assay = "RNA",
                                   variable.features.n = n.variable.genes,
                                   vars.to.regress = col_name,
                                   seed.use = random.seed,
                                   verbose = FALSE)
    } else {
    seurat.object <- SCTransform(seurat.object,
                                 assay = "RNA",
                                 variable.features.n = n.variable.genes,
                                 seed.use = random.seed,
                                 verbose = FALSE)
    }

    # check if PCA components exist in Seurat object
    if (is.null(seurat.object@reductions$PCA))

  if (perform.SingleR == TRUE & is.null(ref.data)) {
    print("Using the Human Primary Cell Atlas as a reference")
    bulk_ref <- HumanPrimaryCellAtlasData()

    }
  }
}
