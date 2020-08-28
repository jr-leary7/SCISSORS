#' Prepare scRNA-seq data for reclustering.
#'
#' This function prepares a `Seurat` object for reclustering analysis.
#' The input is a `Seurat` object in any stage of pre-processing, or even a `SingleCellExperiment` object that will be converted to `Seurat` format.
#' The function checks which metadata features (% mitochondrial DNA, cell cycle scores) and assays are present (normalized counts, PCA & t-SNE embeddings),
#' then runs an initial graph-based clustering.
#' @import Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @param seurat.object The object containing the cells you'd like to analyze.
#' @param n.variable.genes The number of variable genes to find at each step. Defaults to 4000.
#' @param n.PC The number of PCs used as input to non-linear dimension reduction and clustering algorithms. Defaults to 30.
#' @param which.dim.reduc (Optional) Which non-linear dimension reduction algorithms should be used? Supports "tsne", "umap", "phate", and "all". Plots will be generated using the t-SNE embedding. Defaults to c("tsne", "umap"), as most users will likely not have `phateR` installed.
#' @param perplexity (Optional) What perplexity value should be used when embedding cells in t-SNE space? Defaults to 30.
#' @param initial.resolution The initial resolution parameter used in the `FindClusters` function. Defaults to 0.3.
#' @param k.val (Optional) The parameter *k* to be used when creating the shared nearest-neighbor graph. Defaults to *k* ~ sqrt(*n*).
#' @param do.plot (Optional) Should the function print a t-SNE plot of your cells to the graphics viewer? Defaults to FALSE.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @export
#' @examples
#' PrepareData(seurat.object, n.variable.genes = 5000, do.plot = TRUE)
#' PrepareData(seurat.object, initial.resolution = .5, k.val = 25, random.seed = 100)
#' @references
#' Stuart *et al* (2019). Comprehensive integration of single-cell data. *Cell*.

PrepareData <- function(seurat.object = NULL,
                        n.variable.genes = 4000,
                        n.PC = 30,
                        which.dim.reduc = c("tsne", "umap"),
                        perplexity = 30,
                        initial.resolution = .3,
                        k.val = NULL,
                        do.plot = FALSE,
                        random.seed = 629) {
  # check inputs & assays present in Seurat object
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object!") }
  # run function
  if (class(seurat.object)[1] == "SingleCellExperiment") {
    # convert object
    print("Converting user-supplied SingleCellExperiment object to Seurat object")
    seurat.object <- as.Seurat(seurat.object, data = NULL)
    # add necessary metadata to calculate % mito & regress it out
    RNA_counts <- colSums(x = seurat.object, slot = "counts")
    feature_counts <- colSums(x = GetAssayData(object = seurat.object, slot = "counts") > 0)
    seurat.object@meta.data$nCount_RNA <- RNA_counts
    seurat.object@meta.data$nFeature_RNA <- feature_counts
    seurat.object[["percent_MT"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-|^mt-")
    # normalize counts and find highly variable genes
    print("Normalizing counts using SCTransform negative-binomial regression")
    seurat.object <- SCTransform(seurat.object,
                                 assay = "RNA",
                                 vars.to.regress = "percent_MT",
                                 variable.features.n = n.variable.genes,
                                 seed.use = random.seed,
                                 verbose = FALSE)
  }
  else if (is.null(seurat.object@assays$SCT) & length(VariableFeatures(seurat.object)) == 0) {
    # check if % mito DNA exists in Seurat object metadata & regress out if so
    if (any(grepl("MT|mt|Mito|mito", colnames(seurat.object@meta.data)))) {
      col_loc <- which(grepl("MT|mt|Mito|mito", colnames(seurat.object@meta.data)))
      col_name <- colnames(seurat.object@meta.data)[col_loc]
      print("Normalizing counts using SCTransform negative-binomial regression")
      seurat.object <- SCTransform(seurat.object,
                                   assay = "RNA",
                                   variable.features.n = n.variable.genes,
                                   vars.to.regress = col_name,
                                   seed.use = random.seed,
                                   verbose = FALSE)
    } else {
      # add % mito and regress out
      seurat.object[["percent_MT"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-|^mt-")  # non-specific to species
      print("Normalizing counts using SCTransform negative-binomial regression")
      seurat.object <- SCTransform(seurat.object,
                                   assay = "RNA",
                                   variable.features.n = n.variable.genes,
                                   seed.use = random.seed,
                                   verbose = FALSE)
    }
  }

  # add cell cycle scores to Seurat object
  seurat.object <- CellCycleScoring(seurat.object,
                                    s.features = cc.genes.updated.2019$s.genes,
                                    g2m.features = cc.genes.updated.2019$g2m.genes,
                                    set.ident = FALSE)

  # check if PCA components exist in Seurat object
  if (is.null(seurat.object@reductions$pca)) {
    print(sprintf("Running PCA with %s principal components using %s highly variable genes",
                  n.PC,
                  n.variable.genes))
    seurat.object <- RunPCA(seurat.object,
                            features = VariableFeatures(seurat.object),
                            npcs = n.PC,
                            verbose = FALSE,
                            seed.use = random.seed)
  }

  # run t-SNE
  if ("tsne" %in% which.dim.reduc) {
    print(sprintf("Running t-SNE on %s principal components with perplexity = %s", n.PC, perplexity))
    seurat.object <- RunTSNE(seurat.object,
                             reduction = "pca",
                             dims = 1:n.PC,
                             dim.embed = 2,
                             seed.use = random.seed,
                             perplexity = perplexity)
  }

  # run UMAP
  if ("umap" %in% which.dim.reduc) {
    print(sprintf("Running UMAP on %s principal components", n.PC))
    seurat.object <- RunUMAP(seurat.object,
                             umap.method = "uwot",
                             dims = 1:n.PC,
                             n.components = 2,
                             reduction = "pca",
                             verbose = FALSE,
                             seed.use = random.seed)
  }

  # run PHATE
  if ("phate" %in% which.dim.reduc) {
    print(sprintf("Running PHATE on %s principal components", n.PC))
    pca_df <- data.frame(Embeddings(seurat.object, reduction = "pca"))
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
    seurat.object@reductions$phate <- phate_obj
  }

  # initial clustering
  # set k if it wasn't user-defined
  if (is.null(k.val)) {k.val <- round(sqrt(ncol(seurat.object)))}
  print(sprintf("Clustering cells in PCA space using k ~ %s & resolution = %s", k.val, initial.resolution))
  seurat.object <- FindNeighbors(seurat.object,
                                 reduction = "pca",
                                 dims = 1:n.PC,
                                 k.param = k.val)
  seurat.object <- FindClusters(seurat.object,
                                resolution = initial.resolution,
                                algorithm = 1,
                                random.seed = 629,
                                verbose = FALSE)
  print(sprintf("Found %s unique clusters", length(unique(seurat.object$seurat_clusters))))

  # plot results, if user desires
  if (do.plot == TRUE) {
    print(DimPlot(seurat.object, reduction = "tsne"))
  }

  return(seurat.object)
}
