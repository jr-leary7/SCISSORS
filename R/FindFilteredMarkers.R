#' This function finds marker genes for a given set of cells, and filters out highly expressed genes from cells in a reference Seurat object.
#'
#' @name FindFilteredMarkers
#' @author Jack Leary
#' @description Used to filter candidate marker genes for a subcluster Seurat object. Prevents markers from being chosen if they are highly expressed in other celltypes.
#' @import magrittr
#' @importFrom Seurat Idents DefaultAssay FindAllMarkers
#' @importFrom dplyr filter mutate group_by summarise across ungroup select pull
#' @importFrom stats quantile
#' @param obj.1
#' @param obj.2
#' @param ident.1
#' @param ident.2
#' @param test.use
#' @param pval.cutoff
#' @param logfc.cutoff
#' @param quantile.cutoff
#' @param extra.cell.filter
#' @seealso \code{\link{FindSpecificMarkers}}
#' @seealso \code{\link[Seurat]{FindAllMarkers}}
#' @export
#' @examples
#' \dontrun{filtered_markers <- FindFilteredMarkers(obj.1 = subclust_obj, obj.2 = full_obj, ident.1 = "label", ident.2 = "seurat_clusters")}


FindFilteredMarkers <- function(obj.1 = NULL,
                                obj.2 = NULL,
                                ident.1 = "seurat_clusters",
                                ident.2 = "seurat_clusters",
                                test.use = "wilcox",
                                pval.cutoff = 0.05,
                                logfc.cutoff = 0.25,
                                quantile.cutoff = 0.9,
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
                                    logfc.threshold = logfc.cutoff,
                                    test.use = test.use,
                                    only.pos = TRUE,
                                    verbose = FALSE,
                                    random.seed = 312) %>%
             dplyr::filter(p_val_adj < pval.cutoff)
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
                     dplyr::filter(V1 > stats::quantile(V1, quantile.cutoff)) %>%
                     dplyr::mutate(gene = rownames(.)) %>%
                     dplyr::pull(gene)
    high_exp_genes <- c(high_exp_genes, top_exp_genes)
  }
  high_exp_genes <- unique(high_exp_genes)
  markers %<>% dplyr::filter(!gene %in% high_exp_genes)
  return(markers)
}
