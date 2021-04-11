#' A function to perform nonlinear dimension reduction using PHATE.
#'
#' This function wraps the PHATE dimension reduction algorithm in the typical Seurat syntax.
#' @import phateR
#' @import Seurat
#' @param object The Seurat object you'd like to run PHATE on. Defaults to NULL.
#' @param n.components How dimensions should the data be reduced to? Defaults to 2.
#' @param n.PC How many PCs should be used as input to the algorithm? Defaults to NULL.
#' @param mds.method The solver used for MDS. Defaults to SMACOF, but SGD can be used to speed up the algorithm.
#' @param dist.metric The distance metric to use for KNN and MDS. Defaults to the cosine distance.
#' @param random.seed The random seed used to control stochasticity. Defaults to 629.
#' @export
#' @examples
#' RunPhate(object = pbmc3k, n.components = 2, n.PC = 10)

RunPHATE <- function(object = NULL,
                     n.components = 2,
                     n.PC = NULL,
                     mds.method = "smacof",
                     dist.metric = "cosine",
                     random.seed = 629) {
  # check inputs
  if (any(sapply(c(object, n.PC), is.null))) stop("Please provide a Seurat object and a number of PCs to use")

  # run PHATE
  pca_df <- data.frame(Embeddings(object, reduction = "pca"))[, 1:n.PC]
  phate_res <- phate(pca_df,
                     ndim = n.components,
                     mds.solver = mds.method,
                     knn.dist.method = dist.metric,
                     mds.dist.method = dist.metric,
                     npca = NULL,
                     seed = random.seed,
                     verbose = FALSE)
  phate_obj <- CreateDimReducObject(embeddings = phate_res$embedding,
                                    assay = "SCT",
                                    key = "PHATE_",
                                    global = TRUE)
  object@reductions$phate <- phate_obj

  return(object)
}
