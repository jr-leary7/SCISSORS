
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `SCISSORS`

<!-- badges: start -->
<!-- badges: end -->

This package enables researchers to efficiently, accurately, and
reproducibly uncover rare celltypes and other small subgroups using an
iterative reclustering optimization routine built as an extension to a
typical `Seurat` workflow. Tools are included to preprocess data,
identify potential targets for subgroup analysis, optimize reclustering,
and annotate results. A preprint detailing the method & some biological
results is available [on
BioRxiv](https://doi.org/10.1101/2021.10.29.466448).

## Installation

You can install the most recent version of `SCISSORS` with:

``` r
remotes::install_github("jr-leary7/SCISSORS")
```

# Usage

## Libraries

First we’ll need to load our package as well as a couple dependencies.

``` r
library(dplyr)
library(Seurat)
library(ggplot2)
library(SCISSORS)
```

## Preprocessing

Following scRNA-seq tool development tradition, we’ll use the [10X
Genomics 3,000 PBMCs
dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)
as an example of our workflow. The raw reads can be downloaded from
10X’s site, and processed counts are available in the `SeuratData`
package.

``` r
seu_pbmc <- SeuratData::LoadData("pbmc3k") %>% 
            PrepareData(n.HVG = 4000, 
                        n.PC = 15, 
                        which.dim.reduc = "umap", 
                        initial.resolution = 0.4, 
                        use.parallel = FALSE, 
                        random.seed = 629)
#> [1] "Running UMAP on 15 principal components"
#> [1] "Found 6 unique clusters"
```

Visualizing the broad clusters on our UMAP embedding, we can see a few
clusters that seem heterogeneous enough to be reclustered.

``` r
DimPlot(seu_pbmc) + 
  labs(x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Louvain\nCluster") + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank())
```

<img src="man/figures/README-umap-broad-cluster-1.png" width="100%" />

We estimate silhouette scores for each cell, then visualize the
distribution for each broad cluster. Clusters 0, 1, & 3 have decently
lower median scores than the other clusters.

``` r
sil_score_df <- ComputeSilhouetteScores(seu_pbmc, avg = FALSE)
ggplot(sil_score_df, aes(x = Cluster, y = Score, fill = Cluster)) + 
  geom_violin(draw_quantiles = 0.5, 
              scale = "width", 
              size = 1, 
              show.legend = FALSE) + 
  labs(x = "Louvain Cluster", 
       y = "Silhouette Score") + 
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_line(color = "grey80"))
```

<img src="man/figures/README-silhouette-dist-1.png" width="100%" />

Looking at the expression of canonical marker genes allows us to assign
a most-likely celltype identity to each broad cluster. This biological
knowledge aids in the reclustering process i.e., we know that it’s
probably best to recluster clusters 0 & 3 simultaneously, as they’re
both composed of T cells. It also gives us some idea of which / how many
cell subtypes to expect in each cluster, which can improve our
confidence in the final results. For example, if we saw 10 subclusters
in the T cell group it wouldn’t make much sense, as it’s not likely that
so many T cell subtypes exist in such a small dataset.

``` r
markers <- c("CD3G", "IL7R",    # CD4+ T
             "LYZ", "CD14",     # CD14+ monocyte
             "MS4A1", "CD79A",  # B 
             "CD8A",            # CD8+ T
             "RHOC",            # CD16+ monocyte
             "NKG7")            # NK
VlnPlot(seu_pbmc, 
        features = markers, 
        stack = TRUE, 
        flip = TRUE, 
        fill.by = "ident") + 
  labs(y = "Expression", 
       fill = "Louvain\nCluster") + 
  theme(axis.title.x = element_blank())
```

<img src="man/figures/README-broad-markers-1.png" width="100%" />

We add the broad celltype labels to our clusters, then visualize the
results.

``` r
seu_pbmc@meta.data <- mutate(seu_pbmc@meta.data, 
                             broad_celltype = case_when(seurat_clusters == "0" ~ "CD4+ T", 
                                                        seurat_clusters == "1" ~ "CD14+ Monocyte", 
                                                        seurat_clusters == "2" ~ "B", 
                                                        seurat_clusters == "3" ~ "CD8+ T", 
                                                        seurat_clusters == "4" ~ "CD16+ Monocyte", 
                                                        seurat_clusters == "5" ~ "NK", 
                                                        TRUE ~ NA_character_), 
                             broad_celltype = factor(broad_celltype, levels = c("CD4+ T", 
                                                                                "CD8+ T", 
                                                                                "NK", "B", 
                                                                                "CD14+ Monocyte", 
                                                                                "CD16+ Monocyte")))
DimPlot(seu_pbmc, group.by = "broad_celltype") + 
  labs(x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Broad Celltype") + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_blank())
```

<img src="man/figures/README-broad-celltypes-1.png" width="100%" />

## Reclustering

Now that we have an idea of which clusters are poorly fit, as well as a
general idea of which celltypes are which, we’re ready to do some
subclustering. We’ll start with the T cells; since we know that CD4+ and
CD8+ T cells are fairly similar, we’ll recluster them together.

``` r
t_reclust <- ReclusterCells(seu_pbmc, 
                            which.clust = c(0, 3), 
                            merge.clusters = TRUE, 
                            k.vals = c(30, 40, 50), 
                            resolution.vals = c(.2, .3, .4), 
                            n.HVG = 4000, 
                            n.PC = 15, 
                            use.parallel = FALSE, 
                            redo.embedding = TRUE, 
                            random.seed = 312)
#> [1] "Reclustering cells in clusters 0, 3 using k = 50 & resolution = 0.3; S = 0.263"
DimPlot(t_reclust) + 
  labs(x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Louvain\nSubcluster") + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank())
```

<img src="man/figures/README-t-reclust-1.png" width="100%" />

Some canonical T-cell subtype markers help us to differentiate between
the subclusters. Cluster 0 houses Naive CD4+ T-cells, and Memory CD4+
T-cells are located in cluster 1. The CD8+ T-cells are split into
effector & memory subpopulations in clusters 2 and 3, respectively.

``` r
tcell_markers <- c("CCR7", "S100A4", "IL7R", "CD44",  # CD4+ T
                   "CD8B", "GZMB", "FGFBP2", "CD8A",  # effector CD8+ T
                   "GZMK", "NCR3", "KLRB1", "AQP3")   # memory CD8+ T
VlnPlot(t_reclust, 
        features = tcell_markers, 
        stack = TRUE, 
        fill.by = "ident", 
        flip = TRUE) + 
  labs(y = "Expression", 
       fill = "Louvain\nSubcluster") + 
  theme(axis.title.x = element_blank())
```

<img src="man/figures/README-reclust-markers-1.png" width="100%" />

After integrating the subcluster labels back into the full dataset with
`IntegrateSubclusters()`, we add new finer-resolution celltype labels to
the cells & visualize the results on the original UMAP embedding.

``` r
seu_pbmc <- IntegrateSubclusters(seu_pbmc, reclust.results = t_reclust)
seu_pbmc@meta.data <- mutate(seu_pbmc@meta.data, 
                             fine_celltype = case_when(seurat_clusters == "1" ~ "CD14+ Monocyte", 
                                                       seurat_clusters == "2" ~ "B", 
                                                       seurat_clusters == "3" ~ "Memory CD8+ T", 
                                                       seurat_clusters == "4" ~ "CD16+ Monocyte", 
                                                       seurat_clusters == "5" ~ "NK", 
                                                       seurat_clusters == "6" ~ "Naive CD4+ T", 
                                                       seurat_clusters == "7" ~ "Memory CD4+ T", 
                                                       seurat_clusters == "8" ~ "Effector CD8+ T", 
                                                       TRUE ~ NA_character_))
DimPlot(seu_pbmc, group.by = "fine_celltype") + 
  labs(x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Fine Celltype") + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_blank())
```

<img src="man/figures/README-integrate-reclust-1.png" width="100%" />

Lastly, we can use the `FindSpecificMarkers()` function to identify
potential marker genes that are specific (in the statistical sense) to
each cluster, which we accomplish by filtering the set of all DE genes
against a list of genes that are highly expressed in other clusters.
This can be helpful in putting together gene sets for pathway analyses,
celltype classifiers, etc., though it does have the drawback of
filtering out some canonical markers that are highly expressed in more
than one celltype - for example **FCGR3A**, which is highly expressed in
both CD16+ monocytes and NK cells, is removed in these DE results. We
pull the top 3 marker genes by mean log2FC for each celltype, then plot
their expression. Note: the results from this function are pre-filtered
to only include genes with adjusted *p*-values lower than some target
threshold; the default is 0.05.

``` r
de_specific <- FindSpecificMarkers(seu_pbmc, 
                                   ident.use = "fine_celltype", 
                                   perc.cutoff = 0.95)
de_specific_top3 <- de_specific %>% 
                    mutate(cluster = as.character(cluster)) %>% 
                    arrange(cluster, desc(avg_log2FC)) %>% 
                    with_groups(cluster, 
                                slice_head, 
                                n = 3)
DotPlot(seu_pbmc, 
        assay = "RNA", 
        features = unique(de_specific_top3$gene), 
        dot.scale = 8, 
        group.by = "fine_celltype") + 
  coord_flip() + 
  theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95), 
        axis.title = element_blank())
```

<img src="man/figures/README-spec-markers-1.png" width="100%" />

# Conclusions & Best Practices

There’s certainly more that can be done with this dataset - the monocyte
clusters have subgroups as well (DCs, intermediate monocytes, even some
platelets). A further tutorial exploring those cells can be found
[here](https://jr-leary7.github.io/quarto-site/tutorials/SCISSORS_Reclustering.html).
This introduction should provide a solid start though, and has hopefully
shown that reclustering analysis can be efficient and easy, while still
providing results we can be confident in. In general, it’s best if
reclustering analyses are semi-supervised; you want to have some idea of
what’s going on / what could be possible biologically, but it’s
important to not have a predefined conclusion that you’re looking for as
well. Using biological background information as well as cluster fit
heuristics to choose reclustering targets is, empirically, a good way to
avoid false positives and ensure that your reclustering results are
biologically meaningful & reproducible. For more complex analyses using
`SCISSORS` see our manuscript, and don’t hesitate to reach out if you’d
like help with using the package.

# Contact Information

This package is developed and maintained by Jack Leary. Feel free to ask
for help via [opening an
issue](https://github.com/jr-leary7/SCISSORS/issues) or via email
(<jrleary@live.unc.edu>) if more detailed assistance is needed.

# References

1.  Leary, J. *et al*. [Sub-cluster identification through
    semi-supervised optimization of rare cell silhouettes (SCISSORS) in
    single-cell sequencing](https://doi.org/10.1101/2021.10.29.466448).
    *BioRxiv* (2021).

2.  10X Genomics. [3k PBMCs from a healthy donor: single cell gene
    expression dataset by Cell Ranger
    v1.1.0](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k).
    *10X Genomics Documentation* (2016).

3.  Rousseeuw, P. [Silhouettes: A graphical aid to the interpretation
    and validation of cluster
    analysis](https://doi.org/10.1016/0377-0427%2887%2990125-7).
    *Journal of Computational and Applied Mathematics* (1987).
