# YehLabClust
This package allows the reclustering of single cell clusters based on the identification of the highly variable genes that are unique to each cluster. Here we'll briefly describe the workflow on a dataset taken from Baron *et al* (2016), which is available as part of the `scRNAseq` package. 

# Data
First, we load the data.
```{r}
library(Seurat)
library(scRNAseq)
```
# Analysis
Next, we use the `ClusterCells` function to automatically convert the data, which is contained in a `SingleCellExperiment` object, to the `Seurat` format. This functions also calculates the percentage of mitochondrial DNA for each cell, then uses `SCTransform` to select highly variable genes, normalize and scale the counts, regress out the effect of the percentage of mitochondrial DNA, run PCA & t-SNE, cluster the cells, and finally run the re-clustering agorithm on each identified cluster. 
```{r}
baron <- scRNAseq::BaronPancreasData()
baron <- ClusterCells(seurat.object = baron, 
                      n.variable.genes = 4000, 
                      initial.resolution = .5, 
                      random.seed = 629)
```

# Results
Once I generate figures, the results will go here! 
