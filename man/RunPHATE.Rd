% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunPHATE.R
\name{RunPHATE}
\alias{RunPHATE}
\title{Perform nonlinear dimension reduction using PHATE.}
\usage{
RunPHATE(
  object = NULL,
  n.components = 2,
  n.PC = NULL,
  mds.method = "smacof",
  dist.metric = "cosine",
  random.seed = 312
)
}
\arguments{
\item{object}{The Seurat object you'd like to run PHATE on. Defaults to NULL.}

\item{n.components}{How dimensions should the data be reduced to? Defaults to \eqn{p = 2}.}

\item{n.PC}{How many PCs should be used as input to the algorithm? Defaults to NULL.}

\item{mds.method}{The solver used for MDS. Defaults to SMACOF, but SGD can be used to speed up the algorithm.}

\item{dist.metric}{The distance metric to use for KNN and MDS. Defaults to the cosine distance.}

\item{random.seed}{The random seed used to control stochasticity. Defaults to 312.}
}
\value{
A \code{Seurat} object with an added dimension reduction object for the \emph{p}-dimensional PHATE embedding.
}
\description{
This function wraps the PHATE dimension reduction algorithm in the typical Seurat syntax.
}
\examples{
\dontrun{
RunPhate(object = pbmc,
         n.components = 2,
         n.PC = 10)
RunPhate(object = pbmc,
         n.components = 5,
         n.PC = 30,
         dist.metric = "euclidean")
}
}
\seealso{
\code{\link[phateR]{phate}}
}
\author{
Jack Leary
}
