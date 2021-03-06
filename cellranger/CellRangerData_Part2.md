---
title: "CellRanger Data - Part 2"
author: "Craig Duncan"
date: "11-14 March 2021; 6 April 2021"
output: html_notebook
---

# Turning raw data into a count matrix

Here is an example of a dataset of blood cell data:

[SatijaTute](https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html)

The raw data for this investigation of human genome is available here (gzip file).  

[Sample](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

When these are unpacked, you will find a hg19 folder with these files:
 - barcodes.tsv
 - genes.tsv
 - matrix.mtx

 It turns out that this is a very common set of files for single cell RNA-seq datasets, as discussed here:
 
 [hbc](https://hbctraining.github.io/scRNA-seq/lessons/readMM_loadData.html)

(those files, or similar, are produced by the 'CellRanger program used by 10X genomics with its Chromium data')

(that also appears to be using very similar barcodes and gene lists to what the Seurat tutorial uses).  In fact the hbctraining github site is an *excellent* introduction to single cell analysis.  There are some links to other good articles, like this: https://hbctraining.github.io/scRNA-seq/lessons/readMM_loadData.html

## (cellular) Barcodes

This tab-separated-values file is a list of barcodes that determine which *cell* the reads originated from.

barcodes.tsv

This is a file full of barcodes.  Each one is 14 letters long, followed by '-1'.

They look like this:

AAACATACAACCAC-1

The file is arranged alphabetically (the last few codes start with TTT)

To quickly count how many lines (in this case, entries) enter this at console

```
echo $(cat barcodes.tsv|wc -l)
```

Ans: 2700

In the setup, CellRanger seems to need a 'fixed list of known-good barcode sequences.', from which it will show the unfiltered list, or the filtered list of only those that are detected.

Ref:
[10x_2](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)

## genes file

Console output:

```
$ echo $(cat genes.tsv|wc -l)
32738
```

This file is of this form:

ENSG00000243485	MIR1302-10
ENSG00000237613	FAM138A
ENSG00000186092	OR4F5
ENSG00000238009	RP11-34P13.7
ENSG00000239945	RP11-34P13.8
ENSG00000237683	AL627309.1
ENSG00000239906	RP11-34P13.14
ENSG00000241599	RP11-34P13.9
ENSG00000228463	AP006222.2
ENSG00000237094	RP4-669L17.10

The explanation given for the genes in a  'features.tsv' file at 10X is:

"For Gene Expression data, the ID corresponds to gene_id in the annotation field of the reference GTF. Similarly, the name corresponds to gene_name in the annotation field of the reference GTF. If no gene_name field is present in the reference GTF, gene name is equivalent to gene ID."

[10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)

You can see the importance of the presence of suitable GTF reference files if you are attempting to perform the alignment and count with your own custom references files.

Genes.tsv appears to be the same, with only gene data appearing (as a result of options set in the 'count' process?)

## matrix (unique molecular identifier or UMI count matrix)

This determines which *transcript molecule* the read originated from.  It will allow identification of PCR duplicates.

The data may contain duplicated read sequences but you know for sure it is a PCR duplicate if the UMIs are the same (since UMIs are unique per molecule *before* PCR takes place.)

The matrix file is the largest -> 28MB

Console: 

```
$ echo $(cat matrix.mtx|wc -l)
2286887
```

"The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column)."

This file is of this form:

```
"%%MatrixMarket matrix coordinate real general
%
32738 2700 2286884
32709 1 4
32707 1 1
32706 1 10
32704 1 1
32703 1 5
32702 1 6
32700 1 10
32699 1 25
32698 1 3
32697 1 8"
```

The first line shows the number of genes (32738), the barcodes(2700), and then the lines in the matrix (3 lines are used for the headings and these totals, so there are 2286884 entries, 3 less than total number of lines)

This file is, however, different to the sparce matrix described in hbctraining.github.io/scRNA-seq above.
(Harvard Chan Bio group: very helpful)

