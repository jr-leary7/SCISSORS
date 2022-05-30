% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindFilteredMarkers.R
\name{FindFilteredMarkers}
\alias{FindFilteredMarkers}
\title{This function finds marker genes for a given set of cells, and filters out highly expressed genes from cells in a reference Seurat object.}
\usage{
FindFilteredMarkers(
  obj.1 = NULL,
  obj.2 = NULL,
  ident.1 = "seurat_clusters",
  ident.2 = "seurat_clusters",
  de.method = "wilcox",
  fdr.cutoff = 0.05,
  log2fc.cutoff = 0.25,
  perc.cutoff = 0.9,
  extra.cell.filter = NULL
)
}
\arguments{
\item{obj.1}{The object containing clusters that you want marker genes for. Defaults to NULL.}

\item{obj.2}{The object containing reference data against which we'll filter the markers in \code{obj.1}. Defaults to NULL.}

\item{ident.1}{The group identifier for \code{obj.1}. Defaults to "seurat_clusters".}

\item{ident.2}{The group identifier for \code{obj.2}. Defaults to "seurat_clusters".}

\item{de.method}{The differential expression method used in \code{\link[Seurat]{FindAllMarkers}}. Defaults to "wilcox".}

\item{fdr.cutoff}{The cutoff used to remove DE genes with non-significant adjusted \emph{p}-values. Defaults to 0.05.}

\item{log2fc.cutoff}{The log2FC cutoff used, in part, to determine whether a gene is differentially expressed. Defaults to 0.25.}

\item{perc.cutoff}{The percentile cutoff used to find highly expressed genes in other cluster. Defaults to 0.9.}

\item{extra.cell.filter}{An optional list of extra cells to filter out of \code{obj.2} other than the cells in \code{obj.1}. Defaults to NULL.}
}
\description{
Used to filter candidate marker genes for a subcluster Seurat object. Prevents markers from being chosen if they are highly expressed in other celltypes.
}
\examples{
\dontrun{filtered_markers <- FindFilteredMarkers(obj.1 = subclust_obj, obj.2 = full_obj, ident.1 = "label", ident.2 = "seurat_clusters")}
}
\seealso{
\code{\link{FindSpecificMarkers}}

\code{\link[Seurat]{FindAllMarkers}}
}
\author{
Jack Leary
}