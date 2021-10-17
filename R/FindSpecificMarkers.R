#' This functions finds specific marker genes for a all clusters in a \code{Seurat} object.
#'
#' @name FindSpecificMarkers
#' @description This function finds marker genes for all clusters, and then filters those markers on a per-cluster basis against the most highly expressed genes in other clusters.
#' @importFrom dplyr mutate group_by summarise across filter pull bind_rows `%>%`
#' @importFrom Seurat FindAllMarkers
#' @param seurat.object The \code{Seurat} object containing clusters for which you'd like marker genes identified. Defaults to NULL.
#' @param ident.use The cell identity to group by. Defaults to "seurat_clusters".
#' @param de.method The differential expression method used in \code{FindAllMarkers()}. Defaults to "wilcox".
#' @param perc.cutoff The percentile cutoff used to find highly expressed genes in other cluster. Defaults to 0.9.
#' @param log2fc.cutoff The log2FC cutoff used, in part, to determine whether a gene is differentially expressed. Defaults to 0.25.
#' @param fdr.cutoff The cutoff used to remove DE genes with non-significant adjusted \emph{p}-values. Defaults to 0.05.
#' @export
#' @examples
#' FindSpecificMarkers(seurat_object, method = "wilcox")

FindSpecificMarkers <- function(seurat.object = NULL,
                                ident.use = "seurat_clusters",
                                de.method = "wilcox",
                                perc.cutoff = 0.9,
                                log2fc.cutoff = 0.25,
                                fdr.cutoff = 0.05) {
  # check inputs
  if (is.null(seurat.object)) { stop("You forgot to provide a Seurat object!") }
  # get highly expressed genes for each cluster
  gene_means_by_clust <- t(seurat.object@assays$SCT@data) %>%
                         as.data.frame() %>%
                         dplyr::mutate(cell_ident = seurat.objec[[ident.use]]) %>%
                         dplyr::group_by(cell_ident) %>%
                         dplyr::summarise(across(where(is.numeric), mean))
  # find list of genes w/ mean expression above 90th percentile of expression in each celltype
  high_exp_genes <- c()
  cluster_labels <- c()
  for (i in gene_means_by_clust$cell_ident) {
    top_exp_genes <- gene_means_by_clust %>%
                     dplyr::filter(cell_ident == i) %>%
                     dplyr::select(-cell_ident) %>%
                     t() %>%
                     as.data.frame() %>%
                     dplyr::filter(V1 > quantile(V1, perc.cutoff)) %>%
                     dplyr::mutate(gene = rownames(.)) %>%
                     dplyr::pull(gene)
    high_exp_genes <- c(high_exp_genes, top_exp_genes)
    cluster_labels <- c(cluster_labels, rep(i, length(top_exp_genes)))
  }
  high_exp_gene_df <- data.frame(High_Exp_Genes = high_exp_genes, Cluster = as.factor(cluster_labels))
  # get marker genes w/ normal method
  marker_genes <- Seurat::FindAllMarkers(seurat.object,
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
    sub_df <- marker_genes %>% dplyr::filter(cluster == i & !(gene %in% outgroup_high_exp_genes))
    specific_marker_genes <- specific_marker_genes %>% dplyr::bind_rows(sub_df)
  }
  return(specific_marker_genes)
}