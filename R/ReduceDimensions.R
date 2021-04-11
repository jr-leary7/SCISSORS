#' A function to run PCA / t-SNE / UMAP / PHATE.
#'
#' This function simplifies the running of various dimension reduction algorithms. It exists mostly to make the body of ReclusterCells() easier to read.
#' @import Seurat
#' @param obj The Seurat object to run dimension reduction algorithms on. Defaults to NULL.
#' @param n.PC How many principal components should be used? Can be an integer or "auto". Defaults to NULL.
#' @param which.algs Which nonlinear dimension algorithms can be used? Should be some combination of "tsne", "umap", and "phate". Defaults to NULL.
#' @param seed The random seed to use. Defaults to NULL.
#' @export
#' @examples
#' ReduceDimensions(pbmc3k, n.PC = 10, which.algs = "umap", seed = 629)
ReduceDimensions <- function(obj = NULL,
                             n.PC = NULL,
                             which.algs = NULL,
                             seed = NULL) {
  # check inputs
  if(any(sapply(c(obj, n.PC, seed), is.null))) stop("Please provide all inputs to ReduceDimensions")
  which.algs <- tolower(which.algs)  # just in case of

  # reduce dimensions
  if (n.PC != "auto") {
    obj <- RunPCA(obj,
                  npcs = n.PC,
                  features = VariableFeatures(obj),
                  seed.use = seed,
                  verbose = FALSE)
  } else {
    obj <- RunPCA(obj,
                  npcs = 30,
                  features = VariableFeatures(obj),
                  seed.use = random.seed,
                  verbose = FALSE)
    n.PC <- ChoosePCs(obj)
  }

  if (!is.null(which.algs)) {
    if ("tsne" %in% which.algs) {
      obj <- RunTSNE(obj,
                     reduction = "pca",
                     dims = 1:n.PC,
                     dim.embed = 2,
                     seed.use = seed)
    }
    if ("umap" %in% which.algs) {
      obj <- RunUMAP(obj,
                     reduction = "pca",
                     dims = 1:n.PC,
                     umap.method = "uwot",
                     n.components = 2,
                     metric = "cosine",
                     seed.use = seed,
                     verbose = FALSE)
    }
    if ("phate" %in% which.algs) {
      obj <- RunPHATE(obj,
                      n.PC = n.PC,
                      random.seed = random.seed)
    }
  }

  return(obj)
}
