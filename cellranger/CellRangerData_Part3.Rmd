---
title: "CellRanger Data - Part 3 (Seurat pipeline)"
author: "Craig Duncan"
date: "11-14 March 2021; 6 April 2021"
output:
  html_document:
    df_print: paged
---

# Seurat

We need at least R version 4 installed to run Seurat.

Having done that, we start to progress with using the above files in the above directory with the pmbc 10X functions.

# Compressed Seurat workflow

The series of commands that can comprise a workflow are here: https://satijalab.org/seurat/articles/essential_commands.html

This workflow is based on the hg19 folder container the 3 matrix files that are the output from CellRanger.

```{R}
# install.packages('Seurat')
library(dplyr)
library(Seurat) 
library(Matrix)
pbmc.counts <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
DimPlot(object = pbmc, reduction = "tsne")
```

The above workflow will crunch the count data, but it doesn't dwell on the preceding data preparation work that is done in the initial sequencing process or for the initial genes-per-cell counting (perhaps using CellRanger software, which internally uses STAR) as contained in the 10X count-data files (the .mtx, barcodes and genes/features files).
