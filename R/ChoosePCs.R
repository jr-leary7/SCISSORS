#' A function to automatically choose the number of principal components.
#'
#' @name ChoosePCs
#' @author Jack Leary
#' @description This function uses the eigenvalues of the principal component matrix to determine the best number of PCs. The default can be chosen automatically, or given by the user. It is intended to be run after \code{\link[Seurat]{RunPCA}}.
#' @importFrom matrixStats rowVars
#' @importFrom Seurat GetAssayData DefaultAssay Stdev
#' @param seurat.obj The object containing our single cell counts and principal component matrix. Defaults to NULL.
#' @param cutoff The cutoff value for cumulative proportion of variance explained. Can be set by the user, or can be determine automatically. Defaults to 15%.
#' @export
#' @examples
#' \dontrun{ChoosePCs(seurat.obj = pbmc, cutoff = .3)}

ChoosePCs <- function(seurat.obj = NULL, cutoff = .15) {
  # check inputs
  if (is.null(seurat.obj)) { stop("Please supply a Seurat object.") }
  # run function
  total_var <- sum(matrixStats::rowVars(Seurat::GetAssayData(seurat.obj, assay = Seurat::DefaultAssay(seurat.obj), slot = "data")))
  eigenvals <- Seurat::Stdev(seurat.obj, reduction = "pca")^2
  prop_var <- eigenvals / total_var
  cum_prop_var <- cumsum(prop_var)
  if (any(cum_prop_var) > cutoff) {
    cutoff_PC <- which.min(cum_prop_var > cutoff)
  } else {
    print(sprintf("Cumulative % of variance explained did not reach cutoff value = %s.", cutoff))
    cutoff_PC <- length(eigenvals)
  }
  return(cutoff_PC)
}
