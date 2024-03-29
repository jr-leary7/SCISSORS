% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotateMarkerGenes.R
\name{AnnotateMarkerGenes}
\alias{AnnotateMarkerGenes}
\title{Annotate differentially expressed genes using \code{biomaRt}.}
\usage{
AnnotateMarkerGenes(
  marker.genes = NULL,
  species = "human",
  desired.annos = NULL
)
}
\arguments{
\item{marker.genes}{The dataframe of marker genes generated by \code{\link{FindSubpopulationMarkers}}.}

\item{species}{The species of the cells being analyzed. Defaults to "human", but also supports "mouse".}

\item{desired.annos}{The vector containing the annotations you'd like to retrieve for each gene.}
}
\description{
This function uses the \code{biomaRt} package to fetch a user-defined list of attributes for a list of dataframes containing differentially expressed genes. Intended to be run directly after \code{\link{FindSubpopulationMarkers}}.
}
\author{
Jack Leary
}
