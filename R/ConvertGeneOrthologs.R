#' Convert a vector of MGI symbols to their HGNC orthologs and vice versa.
#'
#' @name ConvertGeneOrthologs
#' @author Jack Leary
#' @description Converts gene names from human -> mouse and mouse -> human.
#' @importFrom biomaRt useMart getLDS
#' @param gene.vec A vector of genes to convert. Default to NULL.
#' @param species One of "mm" or "hs". Defaults to "mm" (and thus converts to human).
#' @export
#' @examples
#' \dontrun{ConvertGeneOrthologs(gene.vec = mouse_genes)}
#' \dontrun{ConvertGeneOrthologs(gene.vec = human_genes, species = "hs")}

ConvertGeneOrthologs <- function(gene.vec = NULL, species = "mm") {
  # check inputs
  species <- tolower(species)
  # get marts
  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  # convert gene names
  if (species == "mm") {
    genesV2 <- biomaRt::getLDS(attributes = c("mgi_symbol"),
                               filters = "mgi_symbol",
                               values = gene.vec,
                               mart = mouse,
                               attributesL = c("hgnc_symbol"),
                               martL = human,
                               uniqueRows = TRUE)
    new_genes <- unique(genesV2[, 2])
  } else if (species == "hs") {
    genesV2 <- biomaRt::getLDS(attributes = c("hgnc_symbol"),
                               filters = "hgnc_symbol",
                               values = gene.vec,
                               mart = human,
                               attributesL = c("mgi_symbol"),
                               martL = mouse,
                               uniqueRows = TRUE)
    new_genes <- unique(genesV2[, 2])
  }
  return(new_genes)
}
