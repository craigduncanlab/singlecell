---
title: "Seurat_GeneExpression"
author: "Craig Duncan"
date: "06/04/2021"
output: html_document
---

# Gene Expression analysis

By now, we should have a Seurat object with a data matrix (or matrixes) with genes as row names and barcodes as cell/column names.  Also, there is a data.frame with meta-data by cell barcodes as rownames.

The first really significant part of the gene-per-cell-count data analysis (apart from regressing or filtering out the poor quality data) that Seurat intends to help with is creating and marking sets of highly variable genes (relative to the set of genes in the data), for clustering and other visualisations.

All of these approaches are grounded in statistics, but the actual process is a little more involved than looking at 'variation', and it has to work with the thousands (or millions) of cells in the original RNA seq data.  

In this context, the statistical measures and steps include:
1. The average expression (count) for each gene (across all cells) <---bulk RNA seq method
2. The average dispersion for each gene (across all cells). 
3. Creating 'bins' to hold these values.  i.e. groups.  Say 20 of them.
4. Calculating the dispersion of gene expression within each bin. [z-scores, stats on subsets]

The last step is used to identify the most 'highly variable' genes, in that if there are outlier genes that are highly variable even within each bin, they can be identified as highly variable. 

## Reproducible data science?

The Macosko 2015 paper illustrates how 'reproducible' science, extending to the practical task of software-assisted data analysis, was still in development.  In the Macosko paper, the write-up was split into 2 parts - a very brief reference to using Seurat in the main paper (but with no specific reference to functions), and then a more detailed 'extended' paper that described what was done, and the statistical basis, but not how this was integrated with the Seurat program or its workflow.

By the time we have the 2016 Seurat paper, there is another broad reference to this:
```
The method involves calculating the average expression and dispersion for each gene, placing these genes into bins, and then calculating a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko et al.),
```
## Macosko paper - reduced dimensions and 'training set'

Detailed information on the statistical approach is available in Macosko et al. (2015, Cell). 
Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets.  Cell 161, 1202–1214, May 21, 2015 (The lab of Evan Macosko is currently at the Broad Institute: https://macoskolab.com/papers/)

In the Macosko analysis, there were over 44,000 cells in the original RNA seq data.

The method in the Macosko paper involves working out the reduced dimensions on a smaller subset of cells, then projecting the balance of cells onto the two-dimensional map.

## The t-SNE approach (via PC set generation)

The idea is to compare expression levels to every other gene, to form 'classes' of cells (colour-code these) and then plot them in just two dimensions using t-SNE.

The 'stochastic neighbour embedding' is a 2D visualisation technique, based on higher-dimensional points and assigning probabilities to data objects (based on similarity in the higher dimensions, and in lower dimensions).  This technique will often produce clusters in the 2D representation.  To make these clusters meaningful, the parameter setting must be reasonable.  It is a technique used in a variety of settings, not merely bioinformatics.

In the single cell expression context, this may provides a visual mapping (classification scheme) for identifying cell types within the clusters. 

In some contexts, this can be referred to as 'unsupervised learning', in that optimisation of the parameters in statistical models (like the t-distribution) are used to minimise the divergence between points at a stage in the analysis.   

Crude application of the t-SNE approach will not work if there is too much noise in the data, or there is insufficient information to produce useful models.  In the context of single cell analysis, one of the choices that has been made by researchers (e.g. Satija lab) is to work with a training set of cells in the expression matrix first, which is based on richer information (i.e. higher gene expression counts), then to use this to model to the t-SNE stage, then to proceed with further analysis by projecting all cells into this t-SNE map.

The t-SNE output is the result of a process involving PC sets on a smaller number of cells, then reducing the dimensions, then projecting the balance of the cells into this t-SNE map.

## PCA in general

""

[ref](https://medium.com/analytics-vidhya/eigenvectors-and-eigenvalues-and-there-use-in-principal-component-analysis-machine-learning-1f97fdbdb303)

There is a trade off in accuracy, for simplicity but we try to ensure that it is not enough to worry about.

One simple example is the old school task of walking in a straight line, and then measuring distance to notable objects at a perpendicular distance from your origin line.

The goal of PCA in 2D is to identify the origin line which optimises the data (the least walking distance).  It is a little like linear regression: minimising the mean of point distances from the chosen line. (or using the points to derive this line).  

These can variously be described as mathematical models, statistical modelling. When the parameters or end goal is unknown but can be found by iteration, it is called machine learning (of the origin line).

[The relationship with eigenvectors is that the origin line represents a direction that can simply 'scale' a lot of the points, as an eigenvector does, and only has to deal with a minimal number of points off that line.].  Principal components are, by nature, direction-finding lines for the data, because they aim to explain as much of the data without variance from the line, in a lower dimension than the original data

PCA in general does this for whatever dimension we start with.   For 3D data, we can explain it in 2D principal components.  This is the 'PC subspace representation' of the original data.

## What preparation to data is done before PC?

Log-normalisation is one common method.

Simple, library-size normalisation is intended to standarise the data with respect to the size of libraries that affect the resultant count data (possibly influenced by PCR amplification).   The normalisation adjustments are intended to state each cell's count on a standard scale, with respect to the total counts in the library.  This 'normalises' libraries of different size.  

You might think that this is not necessary if you have only 1 library, but what this also does is standardise the way you present your data, so that it is easily comparable to other studies that use the same technique.  If all your data is normalised so that it has a mean size factor (across all cells) of '1', then you know this is proportionately comparable to another library of cell data that has also been normalised in the same way.

Library size normalization does not take into account batch bias (a well known effect), but it is usually for initial applciations to identify clusters and the top markers that define each cluster.  Beyond this, some other data adjustments may be recommended.

See:
1. http://bioconductor.org/books/release/OSCA/normalization.html#motivation
2. http://bioconductor.org/books/release/OSCA/dimensionality-reduction.html#principal-components-analysis

## Which PCs are most important?

Why do we end up with so many 'PC' options for single cell data?  Each cell represents a dimension.  A scatterplot of any given set of [2 genes scatterplot?] provides a method of reducing that to a lower dimension.  In the usual PC analysis tools, each successive PC is less likely to be the main explanatory PC of the data.

See also:

"By definition, the top PCs capture the dominant factors of heterogeneity in the data set. Thus, we can perform dimensionality reduction by restricting downstream analyses to the top PCs. This strategy is simple, highly effective and widely used throughout the data sciences.".  There is also a general assumption that if there is a structured/coordinated biological process, then it affects cells throughout the data set (whereas individual, independent factors are not likely to be reflected in the principal components description of the data).  As a result, noise is more likely to be captured in the later PCs.

http://bioconductor.org/books/release/OSCA/dimensionality-reduction.html#principal-components-analysis

## The initial PC sets

For single cell analysis, Principal Components Analysis (PCA) is performed on a smaller 'training set'.  This is to improve the analysis, but it also reduces the work to be done.  The goal with the Macoska experiment was to find even more highly variable and statistically significant cell types.

The goal was to identify 'primary structures in the data'. By this, it means principal eigenvectors (also “PC subspace representation”).

The initial output of 'principal components' is 'equal to the number of profiled cells', but only a few are significant.  The goal is to then identify which of the possible PC's is explicable for the data, which involves a technique of randomised perturbation.

The PCA method seems to be based on:
1. (Shalek et al.,2013).  
2. (Chung and Storey, 2014), applied to single-cell RNA-seq data (Shalek et al., 2014).
3. joint-null criterion (Leek and Storey, 2011)

This doesn't require Seurat for the final steps, in that they can:
1. Identify highly variable genes using Seurat
2. Scale and centre data along each gene
3. Use the 'prcomp' function in R.
4. Ascertain if there are 'canonical markers' for different cell types along PCs

## basic t-sne

R has a t-sne function, and the supplementary extended section of the Macosko (2015) paper makes mention of it

"t-Distributed Stochastic Neighbor Embedding (tSNE) (van der Maaten and Hinton, 2008), as implemented in the tsne package in R with the “perplexity” parameter set to 30."

Further, the basis for clustering is explained this way:

"Cells with similar expression signatures of genes within our variable set, and therefore similar PC loadings, will likely localize near each other in the embedding, and hence distinct cell types should form two-dimensional point clouds across the tSNE map."

## Clusters are defined by the 2D representation after t-SNE

Using this approach, the step of finding PC's will usually generate two to 3 dozen significant PC's.   These are reduced to just two dimensions use t-SNE (stochastic neighbour embedding).  The clusters that are visible in two dimensions form the basis for the rest of the analysis.

Here's another explanation of the same process in the 2016 'Seurat' paper (page 1217):
```
"We performed principal components analysis on
the 13,155 largest libraries (Figure S5, Table S3), then reduced
the 32 statistically significant PCs (Experimental Procedures)
to two dimensions using t-Distributed Stochastic Neighbor
Embedding (tSNE) (Amir et al., 2013; van der Maaten and Hinton,
2008). We projected the remaining 36,145 cells in the data into
the tSNE analysis. 

We then combined a density clustering
approach with post hoc differential expression analysis to divide
44,808 cells among 39 transcriptionally distinct clusters (Supplemental
Experimental Procedures) ranging from 50 to 29,400
cells in size (Figures 5B and 5C). 

Finally, we organized the 39
cell populations into larger categories (classes) by building a
dendrogram of similarity relationships among the 39 cell populations
(Figure 5D, left).

The cell populations inferred from this analysis were readily
matched to the known retinal cell types, including all five
neuronal cell classes, based on the specific expression of known
markers for these cell types (Figure 5D, right, and Figure S6A).
```
There are supplementary materials that help explain how this analysis was pursued using Seurat.

```
Principal Components and Clustering Analysis of Retina Data
The clustering algorithm for the retinal cell data was implemented and performed
using Seurat, a recently developed R package for single-cell analysis
(Satija et al., 2015). PCA was first performed on a 13,155-cell ‘‘training set’’
of the 49,300-cell dataset, using single-cell libraries in which transcripts from
>900 genes were detected. We found this approach was more effective in
discovering structures corresponding to rare cell types than performing PCA
on the full dataset, which was dominated by numerous, tiny rod photoreceptors
(Supplemental Experimental Procedures). 

Thirty-two statistically significant
PCs were identified using a permutation test and independently
confirmed using a modified resampling procedure (Chung and Storey, 2015).

We projected individual cells within the training set based on their PC scores
onto a single two-dimensional map using t-Distributed Stochastic Neighbor
Embedding (t-SNE) (van der Maaten and Hinton, 2008). 

The remaining
36,145 single-cell libraries (<900 genes detected) were next projected on
this t-SNE map, based on their representation within the PC-subspace of
the training set (Berman et al., 2014; Shekhar et al., 2014). This approach mit
mitigates
the impact of noisy variation in the lower complexity libraries due to
gene dropouts. It was also reliable in the sense that when we withheld from
the t-SNE all cells from a given cluster and then tried to project them, these
withheld cells were not spuriously assigned to another cluster by the projection
(Table S7). Point clouds on the t-SNE map represent candidate cell types; density
clustering (Ester et al., 1996) identified these regions. Differential expression
testing (McDavid et al., 2013) was then used to confirm that clusters
were distinct from each other. Hierarchical clustering based on Euclidean distance
and complete linkage was used to build a tree relating the clusters. We
noted expression of several rod-specific genes, such as Rho and Nrl, in every
cell cluster, an observation that has been made in another retinal cell gene
expression study (Siegert et al., 2012) and likely arises from solubilization
of these high-abundance transcripts during cell suspension preparation.
Additional information regarding retinal cell data analysis can be found in the
Supplemental Experimental Procedures.
```

Some fine-tuning was also evident, by varying the number of cells used for analysis:

```
"We examined how the classification of cells (based on their
patterns of gene expression) evolved as a function of the
numbers of cells in analysis. We used 500, 2,000, or 9,731 cells
from our dataset, and asked how (for example) cells identified as
amacrines in the full dataset clustered in analyses of smaller
numbers of cells (Figure 5F). As the number of cells in the data
increased, distinctions between related clusters become clearer,
stronger, and finer in resolution, with the result that a greater
number of rare amacrine cell sub-populations (each representing
0.1%–0.9% of the cells in the experiment) could ultimately
be distinguished from one another (Figure 5F)."
```

