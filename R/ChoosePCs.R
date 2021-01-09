#' A function to automatically choose the number of principal components.
#'
#' This function uses the eigenvalues of the principal component matrix to determine the best number of PCs. The default can be chosen automatically, or given by the user. It is intended to be run after `RunPCA()`.
#' @import matrixStats
#' @param seurat.object The object containing our single cell counts and principal component matrix. Defaults to NULL.
#' @param cutoff The cutoff value for cumulative proportion of variance explained. Can be set by the user, or can be determine automatically. Defaults to NULL.
#' @export
#' @examples
#' ChoosePCs(seurat.obj = pbmc, cutoff = .15)

ChoosePCs <- function(seurat.obj = NULL, cutoff = NULL) {
  # check inputs
  if (is.null(seurat.obj)) stop("Please supply a Seurat object.")
  if (is.null(cutoff)) cutoff <- .15
  # run function
  total_var <- sum(rowVars(GetAssayData(seurat.obj, assay = "SCT", slot = "scale.data")))
  eigenvals <- Stdev(seurat.obj, reduction = "pca")^2
  prop_var <- eigenvals / total_var
  cum_prop_var <- c()
  for (i in seq(prop_var)) {
    cum_prop_var[i] <- sum(prop_var[1:i])
  }
  if (any(cum_prop_var) > cutoff) {
    cutoff_PC <- min(which(cum_prop_var > cutoff))
  } else {
    print(sprintf("Cumulative % of variance explained did not reach cutoff value = %s.", cutoff))
    cutoff_PC <- length(eigenvals)
  }
  return(cutoff_PC)
}
