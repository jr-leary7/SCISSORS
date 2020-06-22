#' Iteratively recluster cells.
#'
#' This function prepares a `Seurat` object for iterative reclustering. The input is a `Seurat` object in any stage of pre-processing, or even a `SingleCellExperiment` object that will be converted to `Seurat` format. The function checks which metadata features and assays are present (% mitochondrial DNA, normalized counts, PCA & t-SNE embeddings), then runs an initial graph-based clustering. Next, it calls `SingleR`, then iteratively reclusters the cell types in order to reveal subtypes and rare cell populations.
#' @import Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @param seurat.object The object containing the cells you'd like to analyze.
#' @param n.variable.genes The number of variable genes to find at each step. Defaults to 3000.
#' @param initial.resolution The initial resolution parameter used in the `FindClusters` function. Defaults to 0.5.
#' @param run.SingleR Should cell type identification through `SingleR` be run? Defaults to TRUE.
#' @param ref.data (Optional) A user-defined reference dataset to be used with `SingleR`. If NULL, the Human Primary Cell Atlas dataset will be used. Defaults to NULL.
#' @param recluster.res How many times should each cluster be divided? Defaults to 1.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @param ... Optional parameter passed to other functions. Default is NULL.
#' @export
#' @examples
#' ClusterCells(seurat.object = pbmc, n.variable.genes = 5000, initial.resolution = .75, ref.data = cell_ref)
#' ClusterCells(seurat.object = baron_pancreas, initial.resolution = .5, perform.SingleR = TRUE, species = "mouse")
#' @references
#' Stuart *et al* (2019). Comprehensive integration of single-cell data. *Cell*.

ClusterCells <- function(seurat.object = NULL,
                         n.variable.genes = 3000,
                         initial.resolution = .5,
                         random.seed = 629) {
  # check arguments & assays present in Seurat object
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object!") }
  if (class(seurat.object)[1] == "SingleCellExperiment") {
    print("Converting user-supplied SingleCellExperiment object to Seurat object")
    seurat.object <- as.Seurat(seurat.object, data = NULL)
    # add necessary metadata to calculate % mito & regress it out
    RNA_counts <- colSums(x = seurat.object, slot = "counts")
    feature_counts <- colSums(x = GetAssayData(object = seurat.object, slot = "counts") > 0)
    seurat.object@meta.data$nCount_RNA <- RNA_counts
    seurat.object@meta.data$nFeature_RNA <- feature_counts
    seurat.object[["percent_MT"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-|^mt-")
    seurat.object <- SCTransform(seurat.object,
                                 assay = "RNA",
                                 vars.to.regress = "percent_MT",
                                 variable.features.n = n.variable.genes,
                                 seed.use = random.seed,
                                 verbose = FALSE)
  }
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
      # add % mito and regress out
      seurat.object[["percent_MT"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-|^mt-")
      seurat.object <- SCTransform(seurat.object,
                                   assay = "RNA",
                                   variable.features.n = n.variable.genes,
                                   seed.use = random.seed,
                                   verbose = FALSE)
    }
  }

  # check if PCA components exist in Seurat object
  if (is.null(seurat.object@reductions$PCA)) {
    seurat.object <- RunPCA(seurat.object,
                            features = VariableFeatures(seurat.object),
                            npcs = 30,
                            verbose = FALSE,
                            seed.use = random.seed)
  }

  # check if t-SNE components exist in Seurat object
  if (is.null(seurat.object@reductions$tsne)) {
    seurat.object <- RunTSNE(seurat.object,
                             reduction = "pca",
                             dims = 1:30,
                             dim.embed = 2,
                             seed.use = random.seed,
                             perplexity = 30)
  }


  # initial clustering
  seurat.object <- FindNeighbors(seurat.object,
                                 reduction = "pca",
                                 dims = 1:30)
  seurat.object <- FindClusters(seurat.object,
                                resolution = .5,
                                algorithm = 1,
                                random.seed = 629)

  # run SingleR
  if (run.SingleR) {
    seurat.object <- RunSingleR(seurat.object, ref.data = ref.data)
  }

  # recluster each cluster 3 times
  reclust_results <- ReclusterCells(seurat.object, recluster.res = 2)

  # process reclustering results
  ## code goes here !


  return(seurat.object)
}

