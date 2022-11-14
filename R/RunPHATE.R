#' Perform nonlinear dimension reduction using PHATE.
#'
#' @name RunPHATE
#' @author Jack Leary
#' @description This function wraps the PHATE dimension reduction algorithm in the typical Seurat syntax.
#' @importFrom phateR phate
#' @importFrom Seurat Embeddings CreateDimReducObject
#' @param object The Seurat object you'd like to run PHATE on. Defaults to NULL.
#' @param n.components How dimensions should the data be reduced to? Defaults to \eqn{p = 2}.
#' @param n.PC How many PCs should be used as input to the algorithm? Defaults to NULL.
#' @param mds.method The solver used for MDS. Defaults to SMACOF, but SGD can be used to speed up the algorithm.
#' @param dist.metric The distance metric to use for KNN and MDS. Defaults to the cosine distance.
#' @param random.seed The random seed used to control stochasticity. Defaults to 312.
#' @return A \code{Seurat} object with an added dimension reduction object for the \emph{p}-dimensional PHATE embedding.
#' @seealso \code{\link[phateR]{phate}}
#' @export
#' @examples
#' \dontrun{
#' RunPhate(object = pbmc,
#'          n.components = 2,
#'          n.PC = 10)
#' RunPhate(object = pbmc,
#'          n.components = 5,
#'          n.PC = 30,
#'          dist.metric = "euclidean")
#' }

RunPHATE <- function(object = NULL,
                     n.components = 2,
                     n.PC = NULL,
                     mds.method = "smacof",
                     dist.metric = "cosine",
                     random.seed = 312) {
  # check inputs
  if (is.null(object) | is.null(n.PC)) { stop("Please provide a Seurat object and a number of PCs to use to RunPHATE().") }
  # run PHATE on PCs
  pca_df <- data.frame(Seurat::Embeddings(object, reduction = "pca"))[, 1:n.PC]
  phate_res <- phateR::phate(pca_df,
                             ndim = n.components,
                             mds.solver = mds.method,
                             knn.dist.method = dist.metric,
                             mds.dist.method = dist.metric,
                             npca = NULL,
                             seed = random.seed,
                             verbose = FALSE)
  phate_obj <- Seurat::CreateDimReducObject(embeddings = phate_res$embedding,
                                            assay = "SCT",
                                            key = "PHATE_",
                                            global = TRUE)
  object@reductions$phate <- phate_obj

  return(object)
}
