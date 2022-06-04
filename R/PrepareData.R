#' Prepare scRNA-seq data for reclustering.
#'
#' @name PrepareData
#' @author Jack Leary
#' @description This function prepares  single cell data for reclustering analysis. The input is a \code{Seurat} object in any stage of pre-processing, or even a \code{SingleCellExperiment} object that will be converted to \code{Seurat} format. The function checks which metadata features (% mitochondrial DNA, cell cycle scores) and assays are present (normalized counts, PCA & t-SNE embeddings), then runs an initial graph-based clustering.
#' @importFrom SeuratObject as.Seurat colSums
#' @importFrom Seurat GetAssayData VariableFeatures CellCycleScoring PercentageFeatureSet SCTransform NormalizeData FindVariableFeatures ScaleData RunPCA RunTSNE RunUMAP FindNeighbors FindClusters DimPlot
#' @param seurat.object The object containing the cells you'd like to analyze. Defaults to NULL.
#' @param use.sct Should \code{SCTransform} be used for normalization / HVG selection? Defaults to TRUE, otherwise typical log-normalization is used.
#' @param n.HVG The number of highly variable genes to compute. Defaults to 4000.
#' @param regress.mt Should the percentage of mitochondrial DNA be computed and regressed out? Works for mouse / human gene names. Defaults to TRUE.
#' @param regress.cc Should cell cycle scores be computed & regressed out? NOTE: uses human cell cycle genes. Defaults to TRUE.
#' @param n.PC The number of PCs used as input to non-linear dimension reduction and clustering algorithms. Can be chosen by user, or set automatically using \code{\link{ChoosePCs}}. Defaults to "auto".
#' @param var.cutoff (Optional) The proportion of variance explained cutoff to be used when n.PC is set to "auto". Defaults to .15.
#' @param which.dim.reduc (Optional) Which non-linear dimension reduction algorithms should be used? Supports "tsne", "umap", "phate", and "all". Plots will be generated using the t-SNE embedding. Defaults to c("umap"), as most users will likely not have \code{phateR} installed.
#' @param perplexity (Optional) What perplexity value should be used when embedding cells in t-SNE space? Defaults to 30.
#' @param umap.lr (Optional) What learning rate should be used for the UMAP embedding? Defaults to 0.05.
#' @param initial.resolution The initial resolution parameter used in the \code{\link[Seurat]{FindAllMarkers}} function. Defaults to 0.3.
#' @param nn.metric (Optional) The distance metric to be used in computing the SNN graph. Defaults to "cosine".
#' @param k.val (Optional) The nearest-neighbors parameter \emph{k} to be used when creating the shared nearest-neighbor graph. Defaults to \eqn{k \approx \sqrt{n}}.
#' @param do.plot (Optional) The dimension reduction view you'd like plotted. Should be one of "tsne", "umap", "phate", or "pca". Defaults to NULL.
#' @param random.seed The seed used to control stochasticity in several functions. Defaults to 629.
#' @seealso \code{\link{ChoosePCs}}
#' @seealso \code{\link[Seurat]{FindAllMarkers}}
#' @export
#' @examples
#' \dontrun{PrepareData(seurat.object, n.variable.genes = 3000, n.PC = 20, do.plot = TRUE)}
#' \dontrun{PrepareData(seurat.object, initial.resolution = .5, k.val = 25, random.seed = 100)}
#' @references
#' Stuart *et al* (2019). Comprehensive integration of single-cell data. *Cell*.

PrepareData <- function(seurat.object = NULL,
                        use.sct = TRUE,
                        n.HVG = 4000,
                        regress.mt = TRUE,
                        regress.cc = TRUE,
                        n.PC = "auto",
                        var.cutoff = .15,
                        which.dim.reduc = c("umap"),
                        perplexity = 30,
                        umap.lr = 0.05,
                        initial.resolution = .3,
                        nn.metric = "cosine",
                        k.val = NULL,
                        do.plot = NULL,
                        random.seed = 629) {
  # check inputs & assays present in Seurat object
  if (is.null(seurat.object)) { stop("You forgot to supply a Seurat object!") }
  # convert SingleCellExperiment object to Seurat if necessary
  if (class(seurat.object)[1] == "SingleCellExperiment") {
    print("Converting user-supplied SingleCellExperiment object to Seurat object")
    seurat.object <- SeuratObject::as.Seurat(seurat.object, data = NULL)
    # add necessary metadata for normalization
    RNA_counts <- SeuratObject::colSums(x = seurat.object, slot = "counts")
    feature_counts <- colSums(x = Seurat::GetAssayData(object = seurat.object, assay = "RNA", slot = "counts") > 0)
    seurat.object@meta.data$nCount_RNA <- RNA_counts
    seurat.object@meta.data$nFeature_RNA <- feature_counts
  }
  # add cell metadata and normalize
  if (regress.cc) {
    # add cell cycle scores
    seurat.object <- Seurat::CellCycleScoring(seurat.object,
                                              s.features = cc.genes.updated.2019$s.genes,
                                              g2m.features = cc.genes.updated.2019$g2m.genes,
                                              set.ident = FALSE)
    seurat.object$CC_difference <- seurat.object$S.Score - seurat.object$G2M.Score
  }
  if (regress.mt) {
    # add % mitochondrial DNA
    seurat.object[["percent_MT"]] <- Seurat::PercentageFeatureSet(seurat.object, pattern = "^MT-|^mt-")  # works for human & mouse
  }
  # set up variables to regress out
  regression_vars <- c()
  if (regress.cc) {
    regression_vars <- c(regression_vars, "CC_difference")
  }
  if (regress.mt) {
    regression_vars <- c(regression_vars, "percent_MT")
  }
  if (use.sct) {
    # normalize counts
    if (length(regression_vars) > 0) {
      seurat.object <- Seurat::SCTransform(seurat.object,
                                           variable.features.n = n.HVG,
                                           vars.to.regress = regression_vars,
                                           seed.use = random.seed,
                                           verbose = FALSE)
    } else {
      seurat.object <- Seurat::SCTransform(seurat.object,
                                           variable.features.n = n.HVG,
                                           seed.use = random.seed,
                                           verbose = FALSE)
    }
  } else {
    seurat.object <- Seurat::NormalizeData(seurat.object,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           verbose = FALSE)
    seurat.object <- Seurat::FindVariableFeatures(seurat.object,
                                                  selection.method = "vst",
                                                  nfeatures = n.HVG,
                                                  verbose = FALSE)
    if (length(regression_vars) > 0) {
      seurat.object <- Seurat::ScaleData(seurat.object,
                                         vars.to.regress = regression_vars,
                                         model.use = "negbinom",
                                         verbose = FALSE)
    } else {
      seurat.object <- Seurat::ScaleData(seurat.object, verbose = FALSE)
    }
  }

  # dimension reduction - PCA, t-SNE, UMAP, and/or PHATE
  if (is.null(seurat.object@reductions$pca)) {
    if (n.PC != "auto") {
      seurat.object <- Seurat::RunPCA(seurat.object,
                                      features = Seurat::VariableFeatures(seurat.object),
                                      npcs = n.PC,
                                      verbose = FALSE,
                                      seed.use = random.seed)
    } else {
      seurat.object <- Seurat::RunPCA(seurat.object,
                                      features = Seurat::VariableFeatures(seurat.object),
                                      npcs = 100,
                                      verbose = FALSE,
                                      seed.use = random.seed)
      n.PC <- ChoosePCs(seurat.object, cutoff = var.cutoff)
    }
  }
  if ("tsne" %in% which.dim.reduc) {
    print(sprintf("Running t-SNE on %s principal components with perplexity = %s", n.PC, perplexity))
    seurat.object <- Seurat::RunTSNE(seurat.object,
                                     reduction = "pca",
                                     dims = 1:n.PC,
                                     dim.embed = 2,
                                     seed.use = random.seed,
                                     perplexity = perplexity)
  }
  if ("umap" %in% which.dim.reduc) {
    print(sprintf("Running UMAP on %s principal components", n.PC))
    seurat.object <- Seurat::RunUMAP(seurat.object,
                                     umap.method = "uwot",
                                     dims = 1:n.PC,
                                     n.components = 2,
                                     learning.rate = umap.lr,
                                     reduction = "pca",
                                     verbose = FALSE,
                                     seed.use = random.seed)
  }
  if ("phate" %in% which.dim.reduc) {
    seurat.object <- RunPHATE(seurat.object,
                              n.components = 2,
                              n.PC = n.PC,
                              random.seed = random.seed)
  }
  # initial clustering
  if (is.null(k.val)) k.val <- round(sqrt(ncol(seurat.object)))  # set k if not defined
  seurat.object <- Seurat::FindNeighbors(seurat.object,
                                         reduction = "pca",
                                         dims = 1:n.PC,
                                         k.param = k.val,
                                         annoy.metric = nn.metric,
                                         nn.method = "annoy",
                                         verbose = FALSE)
  seurat.object <- Seurat::FindClusters(seurat.object,
                                        resolution = initial.resolution,
                                        algorithm = 1,
                                        random.seed = random.seed,
                                        verbose = FALSE)
  print(sprintf("Found %s unique clusters", length(unique(seurat.object$seurat_clusters))))
  # plot results if desired
  if (!is.null(do.plot)) print(Seurat::DimPlot(seurat.object, reduction = do.plot))
  # return prepared object
  return(seurat.object)
}
