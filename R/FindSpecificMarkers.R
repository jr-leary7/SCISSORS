#' Find highly specific marker genes for a all clusters in a \code{Seurat} object.
#'
#' @name FindSpecificMarkers
#' @author Jack Leary
#' @description This function finds marker genes for all clusters, and then filters those markers on a per-cluster basis against the most highly expressed genes in other clusters.
#' @import magrittr
#' @importFrom SeuratObject GetAssayData
#' @importFrom Matrix t
#' @importFrom dplyr mutate group_by summarise across filter pull bind_rows
#' @importFrom Seurat FindAllMarkers
#' @importFrom stats quantile
#' @param seurat.object The \code{Seurat} object containing clusters for which you'd like marker genes identified. Defaults to NULL.
#' @param assay.use The assay to use when testing. Defaults to "RNA".
#' @param slot.use The matrix to use when testing. Defaults to "data", for normalized counts, but "counts" (for the raw counts) can also be used.
#' @param ident.use The cell identity to group by. Defaults to "seurat_clusters".
#' @param de.method The differential expression method used in \code{\link[Seurat]{FindAllMarkers}}. Defaults to "wilcox".
#' @param perc.cutoff The percentile cutoff used to find highly expressed genes in other cluster. Defaults to 0.9.
#' @param log2fc.cutoff The log2FC cutoff used, in part, to determine whether a gene is differentially expressed. Defaults to 0.25.
#' @param fdr.cutoff The cutoff used to remove DE genes with non-significant adjusted \emph{p}-values. Defaults to 0.05.
#' @return A data.frame of markers for each cluster filtered by genes considered highly expressed in other clusters.
#' @seealso \code{\link[Seurat]{FindAllMarkers}}
#' @export
#' @examples
#' \dontrun{FindSpecificMarkers(seurat_object, method = "wilcox", ident.use = "celltype")}
#' \dontrun{FindSpecificMarkers(seurat_object, method = "wilcox", assay.use = "SCT", slot.use = "data")}

FindSpecificMarkers <- function(seurat.object = NULL,
                                assay.use = "RNA",
                                slot.use = "data",
                                ident.use = "seurat_clusters",
                                de.method = "wilcox",
                                perc.cutoff = 0.9,
                                log2fc.cutoff = 0.25,
                                fdr.cutoff = 0.05) {
  # check inputs
  if (is.null(seurat.object)) { stop("You forgot to provide a Seurat object!") }
  if (!ident.use %in% colnames(seurat.object@meta.data)) { stop("ident.use must exist in the object's metadata.") }
  Idents(seurat.object) <- ident.use
  # get highly expressed genes for each cluster
  gene_means_by_clust <- Matrix::t(SeuratObject::GetAssayData(seurat.object, assay = assay.use, slot = slot.use)) %>%
                         as.data.frame() %>%
                         dplyr::mutate(cell_ident = unname(unlist(seurat.object[[ident.use]]))) %>%
                         dplyr::group_by(cell_ident) %>%
                         dplyr::summarise(dplyr::across(where(is.numeric), mean))
  # find list of genes w/ mean expression above 90th percentile of expression in each celltype
  high_exp_genes <- c()
  cluster_labels <- c()
  for (i in gene_means_by_clust$cell_ident) {
    top_exp_genes <- gene_means_by_clust %>%
                     dplyr::filter(cell_ident == i) %>%
                     dplyr::select(-cell_ident) %>%
                     t() %>%
                     as.data.frame() %>%
                     dplyr::filter(V1 > stats::quantile(V1, perc.cutoff)) %>%
                     dplyr::mutate(gene = rownames(.)) %>%
                     dplyr::pull(gene)
    high_exp_genes <- c(high_exp_genes, top_exp_genes)
    cluster_labels <- c(cluster_labels, rep(i, length(top_exp_genes)))
  }
  high_exp_gene_df <- data.frame(High_Exp_Genes = high_exp_genes,
                                 Cluster = as.factor(cluster_labels))
  # get marker genes w/ normal method
  marker_genes <- Seurat::FindAllMarkers(seurat.object,
                                         assay = assay.use,
                                         slot = slot.use,
                                         logfc.threshold = log2fc.cutoff,
                                         test.use = de.method,
                                         only.pos = TRUE,
                                         verbose = FALSE,
                                         random.seed = 312) %>%
                  dplyr::filter(p_val_adj < fdr.cutoff)
  # remove highly expressed genes in other clusters from each cluster's markers
  specific_marker_genes <- NULL
  for (i in unique(marker_genes$cluster)) {
    outgroup_high_exp_genes <- high_exp_gene_df %>%
                               dplyr::filter(Cluster != i) %>%
                               dplyr::pull(High_Exp_Genes) %>%
                               unique()
    sub_df <- dplyr::filter(marker_genes,
                            cluster == i,
                            !(gene %in% outgroup_high_exp_genes))
    specific_marker_genes <- specific_marker_genes %>%
                             dplyr::bind_rows(sub_df)
  }
  return(specific_marker_genes)
}
