#' This function finds marker genes for a given set of cells, and filters out highly expressed genes from cells in a reference Seurat object.
#'
#' @name FindFilteredMarkers
#' @author Jack Leary
#' @description Used to filter candidate marker genes for a subcluster Seurat object. Prevents markers from being chosen if they are highly expressed in other celltypes.
#' @import magrittr
#' @importFrom Seurat Idents DefaultAssay FindAllMarkers
#' @importFrom dplyr filter mutate group_by summarise across ungroup select pull
#' @importFrom stats quantile
#' @param obj.1 The object containing clusters that you want marker genes for. Defaults to NULL.
#' @param obj.2 The object containing reference data against which we'll filter the markers in \code{obj.1}. Defaults to NULL.
#' @param ident.1 The group identifier for \code{obj.1}. Defaults to "seurat_clusters".
#' @param ident.2 The group identifier for \code{obj.2}. Defaults to "seurat_clusters".
#' @param de.method The differential expression method used in \code{\link[Seurat]{FindAllMarkers}}. Defaults to "wilcox".
#' @param fdr.cutoff The cutoff used to remove DE genes with non-significant adjusted \emph{p}-values. Defaults to 0.05.
#' @param log2fc.cutoff The log2FC cutoff used, in part, to determine whether a gene is differentially expressed. Defaults to 0.25.
#' @param perc.cutoff The percentile cutoff used to find highly expressed genes in other cluster. Defaults to 0.9.
#' @param extra.cell.filter An optional list of extra cells to filter out of \code{obj.2} other than the cells in \code{obj.1}. Defaults to NULL.
#' @seealso \code{\link{FindSpecificMarkers}}
#' @seealso \code{\link[Seurat]{FindAllMarkers}}
#' @export
#' @examples
#' \dontrun{filtered_markers <- FindFilteredMarkers(obj.1 = subclust_obj, obj.2 = full_obj, ident.1 = "label", ident.2 = "seurat_clusters")}


FindFilteredMarkers <- function(obj.1 = NULL,
                                obj.2 = NULL,
                                ident.1 = "seurat_clusters",
                                ident.2 = "seurat_clusters",
                                de.method = "wilcox",
                                fdr.cutoff = 0.05,
                                log2fc.cutoff = 0.25,
                                perc.cutoff = 0.9,
                                extra.cell.filter = NULL) {
  # check inputs
  if (is.null(obj.1) || is.null(obj.2)) { stop("You must provide two non-NULL Seurat objects to FindFilteredMarkers().") }
  if (!ident.1 %in% colnames(obj.1@meta.data)) { stop("ident.1 not in obj.1 metadata.") }
  if (!ident.2 %in% colnames(obj.2@meta.data)) { stop("ident.2 not in obj.2 metadata.") }
  Seurat::Idents(obj.1) <- ident.1
  Seurat::Idents(obj.2) <- ident.2
  Seurat::DefaultAssay(obj.1) <- "SCT"
  Seurat::DefaultAssay(obj.2) <- "SCT"
  markers <- Seurat::FindAllMarkers(obj.1,
                                    logfc.threshold = log2fc.cutoff,
                                    test.use = de.method,
                                    only.pos = TRUE,
                                    verbose = FALSE,
                                    random.seed = 312) %>%
             dplyr::filter(p_val_adj < fdr.cutoff)
  # filter obj.1 cells from obj.2
  obj.2 <- obj.2[, !rownames(obj.2@meta.data) %in% rownames(obj.1@meta.data)]
  if (!is.null(extra.cell.filter)) {
    obj.2 <- obj.2[, !rownames(obj.2@meta.data) %in% extra.cell.filter]
  }
  gene_means_by_clust <- t(obj.2@assays$SCT@data) %>%
                         as.data.frame() %>%
                         dplyr::mutate(label_fine = obj.2[[ident.2]][, 1]) %>%
                         dplyr::group_by(label_fine) %>%
                         dplyr::summarise(dplyr::across(where(is.numeric), mean)) %>%
                         dplyr::ungroup()
  high_exp_genes <- c()
  loop_celltypes <- unique(obj.2[[ident.2]][, 1])
  for (i in loop_celltypes) {
    top_exp_genes <- gene_means_by_clust %>%
                     dplyr::filter(label_fine == i) %>%
                     dplyr::select(-label_fine) %>%
                     t() %>%
                     as.data.frame() %>%
                     dplyr::filter(V1 > stats::quantile(V1, perc.cutoff)) %>%
                     dplyr::mutate(gene = rownames(.)) %>%
                     dplyr::pull(gene)
    high_exp_genes <- c(high_exp_genes, top_exp_genes)
  }
  high_exp_genes <- unique(high_exp_genes)
  markers %<>% dplyr::filter(!gene %in% high_exp_genes)
  return(markers)
}
