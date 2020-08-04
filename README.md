# SCISSORS
This package allows the reclustering of single cell clusters based on the identification of the highly variable genes that are unique to each cluster. Here we'll briefly describe the workflow on a dataset taken from Baron *et al* (2016), which is available as part of the `scRNAseq` package. 

# Data
First, we load the data.
```{r}
library(Seurat)
library(scRNAseq)
```
# Analysis
Next, we use the `PrepareData` function to automatically convert the data, which is contained in a `SingleCellExperiment` object, to the `Seurat` format. This function also calculates the percentage of mitochondrial DNA for each cell and uses `SCTransform` to select highly variable genes, normalize and scale the counts, and regress out the effect of the percentage of mitochondrial DNA. It then runs PCA & t-SNE using 30 PCs as an intialization for the t-SNE embedding. Finally, we create a SNN graph using the approximation $k = \sqrt{N}$ and generate a preliminary rough clustering of our cells using Louvain modularity optimization. Basically, the function performs all the necessary pre-processing steps commonly used in `Seurat`. 
```{r}
baron <- scRNAseq::BaronPancreasData()
baron <- PrepareData(seurat.object = baron, 
                     n.variable.genes = 4000, 
                     initial.resolution = .5, 
                     random.seed = 629)
```

# Reclustering
The `ReclusterCells` function performs the actual subpopulation-detection analysis, which is based on tuning the parameters of the Louvain modularity optimization function with the goal of maximizing the mean silhouette score of a given set of parameters. In this case, after investigating the preliminary clustering results, we decide to identify subpopulations in clusters 0, 3, 5 & 6. The `do.plot` argument allows the printing of the optimal reclustering results for eacgh cluster to the graphics viewer in RStudio. 
```{r}
reclust_results <- ReclusterCells(seurat.object = baron, 
                                  n.variable.genes = 4000, 
                                  which.clust = list(0, 3, 5, 6), 
                                  resolution.vals = c(.05, .1, .2, .35), 
                                  do.plot = TRUE)
```

# Identifying Subpopulation Marker Genes
Next, we'd of course like to identify marker genes for our newly discovered subpopulations. This process is implemented in the `FindSubpopulationMarkers` function. It creates a "new" cluster for the identified subpopulation, and uses a one-versus-many approach to determine which genes uniquely identify it. The test used can be user-specified, but defaults to a Wilcox test. Here we find marker genes for each of the identified subpopulations in cluster 6.
```{r}
subpop_markers <- FindSubpopulationMarkers(seurat.object = baron, 
                                           reclust.data = reclust_results, 
                                           which.clust = list(6), 
                                           diff.exp.test = "wilcox", 
                                           log.fc.thresh = 1.75, 
                                           random.seed = 629))
```

The subpopulation markers can be annotated using `biomaRt` or another gene onology package, and can be used to determine, on a transcriptomic level, what about the molecular biology of your subpopulation makes it unique and special. 

# Contact Information
This package is based on the ideas of Xianlu Peng, PhD (Research Asst. Professor, Lineberger Comprehensive Cancer Center, UNC Chapel Hill). The code is written and maintained by Jack Leary (Research Collaborator, Lineberger Comprehensive Cancer Center, UNC Chapel Hill). Jack can be reached on GitHub as well as at jrleary@live.unc.edu for any questions, issues, or bugs. 
