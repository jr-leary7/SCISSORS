# SCISSORS

This package allows the reclustering of single cell clusters based on the identification of the highly variable genes that are unique to each cluster. Here we'll briefly describe the workflow on the classic PBMC3k dataset from 10X Genomics.

# Installation

`SCISSORS` can be installed from this repository like so:

```{r}
devtools::install_github("jr-leary7/SCISSORS")
```

# Getting Started

Here's a basic tutorial on how to get `SCISSORS` up and running (and finding cool cell subgroups for you!).

## Libraries

First, we load the necessary libraries, `Seurat` to contain our analyses and `SeuratData` for our example dataset.

```{r}
library(Seurat)
library(SeuratData)
```

## Preprocessing

Next, we use the `PrepareData()` function to pre-process our data. This function calculates the percentage of mitochondrial DNA for each cell and uses `SCTransform` to select highly variable genes, normalize and scale the counts, and regress out the effect of the percentage of mitochondrial DNA. It then runs PCA, chooses an appropriate number of principal components using a cutoff value for the cumulative proportion of variance explained, and then performs t-SNE using the principal component matrix as an initialization for the embedding. Finally, we create a SNN graph using the approximation *k* \~ sqrt(*N*) and generate a preliminary rough clustering of our cells using Louvain modularity optimization. Basically, the function performs all the necessary pre-processing steps commonly used in `Seurat`.

```{r}
pbmc <- SeuratData::LoadData("pbmc3k")
pbmc <- PrepareData(seurat.object = pbmc3k, 
                    n.variable.genes = 4000, 
                    n.PC = "auto", 
                    cutoff = .15, 
                    which.dim.reduc = c("tsne", "umap"), 
                    initial.resolution = .3, 
                    random.seed = 629)
```

## Reclustering

The `ReclusterCells()` function performs the actual subpopulation-detection analysis, which is based on tuning the parameters of the Louvain modularity optimization function with the goal of maximizing the mean silhouette score of a given set of parameters. In this case, after investigating the preliminary clustering results, we decide to identify subpopulations in clusters 0, 1, and 2. The `do.plot` argument allows the printing of the optimal reclustering results for each cluster to the graphics viewer in RStudio.

```{r}
reclust_results <- ReclusterCells(seurat.object = pbmc3k, 
                                  which.clust = list(0, 1, 2), 
                                  merge.clusters = FALSE, 
                                  n.PC = 10, 
                                  n.variable.genes = 4000, 
                                  which.dim.reduc = "tsne", 
                                  redo.embedding = TRUE, 
                                  resolution.vals = c(.1, .2, .3), 
                                  k.vals = c(10, 20, 30), 
                                  random.seed = 629)
```

## Identifying Subpopulation Marker Genes

Next, we'd of course like to identify marker genes for our newly discovered subpopulations. This process is implemented in the `FindSubpopulationMarkers()` function. It creates a "new" cluster for the identified subpopulation, and uses a one-versus-many approach to determine which genes uniquely identify it. The test used can be user-specified, but defaults to a Wilcox test. Here we find marker genes for each of the identified subpopulations in cluster 0 within the context of the entire sample.

```{r}
subpop_markers <- FindSubpopulationMarkers(seurat.object = pbmc3k, 
                                           reclust.data = reclust_results[[1]], 
                                           which.compare = "all cells", 
                                           diff.exp.test = "wilcox", 
                                           log.fc.thresh = 1.75, 
                                           random.seed = 629))
```

The subpopulation markers can be annotated using `biomaRt` or another gene ontology package, and can be used to determine, on a transcriptomic level, what about the molecular biology of your subpopulation makes it unique. Eventually, I'll finish the `AnnotateMarkerGenes()` function to perform this automatically on the results of `FindSubpopulationMarkers()`, but as of now that is a low-priority fix for the package.

# Contact Information

This package is based on the ideas of Xianlu Peng, PhD (Research Asst. Professor, Lineberger Comprehensive Cancer Center, UNC Chapel Hill). The code is written and maintained by Jack Leary (Research Collaborator, Lineberger Comprehensive Cancer Center, UNC Chapel Hill & University of Florida Department of Biostatistics). Jack can be reached on GitHub as well as at [jrleary\@live.unc.edu](mailto:jrleary@live.unc.edu) for any questions, issues, or bugs.
