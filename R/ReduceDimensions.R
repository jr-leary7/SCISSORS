#' A function to run PCA / t-SNE / UMAP / PHATE.
#'
#' @name ReduceDimensions
#' @author Jack Leary
#' @description This function simplifies the running of various dimension reduction algorithms. It exists mostly to make the body of \code{\link{ReclusterCells}} easier to read.
#' @importFrom Seurat RunPCA VariableFeatures RunTSNE RunUMAP
#' @param obj The Seurat object to run dimension reduction algorithms on. Defaults to NULL.
#' @param n.PC How many principal components should be used? Can be an integer or "auto". Defaults to NULL.
#' @param which.algos Which nonlinear dimension algorithms can be used? Should be some combination of "tsne", "umap", and "phate". Defaults to "umap".
#' @param random.seed The random seed to use. Defaults to 312.
#' @seealso \code{\link{RunPHATE}}
#' @export
#' @examples
#' \dontrun{ReduceDimensions(pbmc3k, n.PC = 10, which.algos = "umap", random.seed = 629)}
ReduceDimensions <- function(obj = NULL,
                             n.PC = NULL,
                             which.algos = "umap",
                             random.seed = 312) {
  # check inputs
  if (is.null(obj) | is.null(n.PC)) { stop("Please provide all inputs to ReduceDimensions().") }
  which.algs <- tolower(which.algs)  # just in case
  # reduce dimensions
  if (n.PC != "auto") {
    obj <- Seurat::RunPCA(obj,
                          npcs = n.PC,
                          features = Seurat::VariableFeatures(obj),
                          seed.use = random.seed,
                          verbose = FALSE)
  } else {
    obj <- Seurat::RunPCA(obj,
                          npcs = 30,
                          features = Seurat::VariableFeatures(obj),
                          seed.use = random.seed,
                          verbose = FALSE)
    n.PC <- ChoosePCs(obj)
  }
  if ("tsne" %in% which.algos) {
    obj <- Seurat::RunTSNE(obj,
                           reduction = "pca",
                           dims = 1:n.PC,
                           dim.embed = 2,
                           seed.use = random.seed)
  }
  if ("umap" %in% which.algos) {
    obj <- Seurat::RunUMAP(obj,
                           reduction = "pca",
                           dims = 1:n.PC,
                           umap.method = "uwot",
                           n.components = 2,
                           metric = "cosine",
                           seed.use = random.seed,
                           verbose = FALSE)
  }
  if ("phate" %in% which.algos) {
    obj <- RunPHATE(obj,
                    n.PC = n.PC,
                    random.seed = random.seed)
  }
  return(obj)
}
