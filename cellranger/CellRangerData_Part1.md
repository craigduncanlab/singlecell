---
title: "CellRanger Data - Part 1"
author: "Craig Duncan"
date: "11-14 March 2021; 6 April 2021"
output: html_notebook
---

# Relevance of CellRanger software to Seurat workflows

This site confirms that Seurat is generally set up to work with CellRanger-specific files (if you have those barcodes etc).  This should be born in mind when starting with some of the Seurat tutorials

[Link](https://adinasarapu.github.io/posts/2019/01/blog-post-sc-ranseq/)

# What is done to obtain the 10X count-data files from say raw base-call sequencing data?

Work is done in the initial sequencing process or the initial genes-per-cell counting (perhaps using CellRanger, which internally uses STAR) that ultimately produces the 10X count-data files (the .mtx, barcodes and genes/features files).

## From FASTA to some kind of count matrix with or without CellRanger?

If you are only given the 'lanes' of FASTQ files, then is that an order of magnitude more difficult to use for a count matrix or not?  And why would that be the case?  What workarounds are necessary?

This site (by Andrew Severin, Iowa) does a pretty good job of explaining how to go from raw data in FASTA files through to raw data in the raw_feature_bc_matrix folder (produced by CellRanger).  (nb features.tsv is renamed to genes.tsv to work with Seurat):

[Link](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0)

There's even a reference to Dave Tang on that page:
[Link](https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/)

On Dave's page he stated that CellRanger is used to generate counts from the FASTQ files produced by the Chromium process.   

"Cell Ranger is the software provided by 10x Genomics to process Chromium single cell 3′ RNA-seq data"

The processing involves a first step:
1. `cellranger mkfastq` runs base call (bcl) files into FASTQ (this is demultiplexing and conversion; it uses a sample_sheet.csv to name the FASTQ files based on libraries and lanes)
2. 'cellranger count' will work with the FASTQ files and produce the 1 directory with 3 files suitable for inputting into Seurat.  A 9 hour process (on a PC?).  The count produces a BAM file as part of the process.

Demultiplexing is only necessary if you push through multiple samples on one lane (in which case you can separate samples before you create the FASTQ files): see https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/processing-raw-scrna-seq-data.html

CellRanger does some quality checking, and if you are getting some poor results, then trimming those poor quality reads or some adapter content might be necessary (but not otherwise).

# The counts process requires reference genomes (generally)

CellRanger works easily with human reference genome, or mouse data:

"Cell Ranger provides pre-built human (hg19, GRCh38), mouse (mm10), and ercc92 reference packages for read alignment and gene expression quantification in cellranger count."

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

If you aren't working with that (e.g. plants), you will need a reference genome and gtd file.  In other words, you need to create your own reference package file:

"To create and use a custom reference package, Cell Ranger requires a reference genome sequence (FASTA file) and gene annotations (GTF file)."  "For the GTF file, genes must be annotated with feature type 'exon' (column 3). -" (as above)

Internally, CellRanger is using STAR, which brings in its own compatibility issues for these FASTA, GTF files.

Some further information is here:
[Link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#gtf)

"The single-nuclei RNA-seq assay captures unspliced pre-mRNA as well as mature mRNA. However, after alignment, cellranger count only counts reads aligned to exons."

# STAR

STAR is a 'splice-aware' alignment program.  STAR aligns reads to a reference genome and produces a BAM file.

Expression level of each gene per cell can be done simplistic counting of reads that overlap with genomic features:
HT-seq or FeatureCounts.  FeatureCounts can count genes, exons, promoter gene bodies etc.

*Also note that STARsolo, integrated with STAR, can do gene quantification for 10X and other droplet-based scRNA-seq.  It can be up to 10 times faster than the CellRanger Count*.

(Here's a 2019 benchmarking paper, discussing BUS, by Melsted et al: 
[Link1](https://www.biorxiv.org/content/10.1101/673285v2.full) and 
[Link2](https://academic.oup.com/bioinformatics/article-abstract/35/21/4472/5487510?redirectedFrom=fulltext)

The barcode approach is designed to gropu together results with similar barcodes (but there may be errors), which is why they use a 'whitelist' of barcodes to verify the barcodes detected, and correct small errors.  Collapsing of duplicate UMIs (assuming they are all PCR duplicates and not actual du0plicates) is naive, but the usual approach because natural variation is considered minimal.

The premise of tools like Kallisto is that precision is not necessary, and results in time-saving.  As Melsted puts it: "Since detailed base-pair alignment is not necessary to generate a count matrix, pseudoalignment to a reference transcriptome8 suffices."

CellRanger requires time and memory to get the job done.  {Specifically, its use of STAR is the main issue.   STAR is even considered to use significant memory for human or mouse (para 5.8.1, https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/processing-raw-scrna-seq-data.html), so something like Wheat will be even more demanding}.  It may not even be feasible for something like Wheat?  How do people working with Wheat overcome these difficulties?

## What Cellranger count does (internally)

"cellranger count takes FASTQ files from cellranger mkfastq and performs alignment, filtering, and UMI counting."

You will need to supply a 'Cell Ranger compatible transcriptome reference'

[Link](https://bioinformatics.uconn.edu/single-cell-rna-sequencing-cell-ranger-2/#)

For quick results, CellRanger not only provides a summary of the count results, but also some further downstream analysis like dimension reduction can be performed using the Cell Ranger R Kit (an R package).

To replicate this without CellRanger it is necessary to know what data is required for the input, and what alternative software is available.

If you can find examples of CellRanger output (as the above site has), then you can practice using this data with downstream programs like Seurat, which whilst CellRanger-supportive, are not exclusive.

## The Cellranger count output

The CellRanger count files are in a particular form.

Two of the output folders produced by CellRanger are matrices files: raw_gene_bc_matrices and filtered_gene_bc_matrices.

You eventually find the ultimate output of '1 directory, 3 files' with the barcodes.tsv, genes.tsv, matrix.mtx.

In effect, the 'features' (or genes) file, and the barcodes files provide the rows and column names, respectively, and the .mtx file contains the count data.  (Downstream steps will integrate these, and use them in this way, but keeping them as sparse matrices if possible).  Converting sparse matrix formats to dense formats (in CSV) may result in big data/disk space problems.  

There is a further detailed explanation about these final steps here:

e.g.
"After these two filtering steps, each observed barcode, UMI, gene combination is recorded as a UMI count in the unfiltered feature-barcode matrix. The number of reads supporting each counted UMI is also recorded in the molecule info file."

[10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview)

This page contains a detailed explanation of the contents of the 3 files (.csv)

[10x_2](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)

The 3 files represent the 'Market Exchange Format (MEX)' which is a method of representing sparse matrices in an abbreviated form [NIST](https://math.nist.gov/MatrixMarket/formats.html)

Both R and Python support MEX format.

## 10x (CellRanger) output to data files, before being imported into R, and then into Seurat

Seurat has a function Read10X(), that can take in the 3 source files that are output (including the .mtx file) and transform all three of them into a matrix object within the R software.  The dGC Matrix is a sparce matrix object in R that has some attributes (i, p, x) to hold the data.  There's a good discussion of the dgC matrix here:
[RBloggers](https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/)

## Why a sparse matrix in R?

After you have the digital information with barcodes and cell sequences, the first step in sorting through the sequencing data is to classify the fragments by cell (i.e. by barcode).  Secondly, you will find that each group of cell barcodes:
(a) will match only a few different genes (a small percentage of the total gene list); and
(b) the sets of genes for each cell will not usually be the same set as other cells.  

The number of *types* of genes matched by a cell is not the same as the count of the absolute number of cells containing those genes. 

The next step, in theory is to prepare a 2D matrix, with the unique genes in the row headings (vertical).  The cells (barcodes) will be column headings (horizontal).   The data will display the count of cells of a given type that match a particular gene.  Assume for now that a 'match' is the presence of a cell with a sequence from a particular gene (it may actually require a minimum number to be significant). 

However, this information is sparse (as demonstrated in the Seurat tutorial).

The biological reason for the sparsity of gene expression data is explained well in the DropSeq paper in Cell (2015). 

The biological explanation for sparsity is this:   As each cell will only match a few genes [from the transcriptome] at any given time the row matches for a given column are sparse.  *Some quality control of this 'count' data is needed; e.g. find only those cells with a certain minimum number of matches, or normalise the data so that it is relative, not absolute.

# Alternatives to CellRanger/STAR for alignment and gene expression counting

STAR and Kallisto-BUStools are two alignment tools for single cell counting.  You need to check if they are splice aware, and whether that affects the outcome.

Alternatives to STAR, for counting, include Subread and Hisat 2.

{Also note that STAR requires reference genome sequences (FASTA) and annotations (GTF) in order to do the alignment with your sequencing data}

## Kallisto

Kallisto with BUStools claims that the memory intensity is much smaller and more consistent across job sizes and there are time savings too:

"Our scRNA-seq workflow is up to 51 times faster than Cell Ranger and up to 4.75 times faster than Alevin. It is also up to 3.5 times faster than STARsolo: a recent version of the STAR aligner adapted for scRNA-seq (Figure 1e, Supplementary Table 2). "

"In benchmarks on the panel described in this paper, kallisto’s running time was comparable to that of the word count (wc) command applied to the FASTQ files,"

There's a kallisto/BUStools workflow here with the use of Seurat.  Not sure about data flows yet.
https://bustools.github.io/BUS_notebooks_R/velocity.html

### Kallisto as a pseudo-aligner

In 10X data, the sequences in 'R1' have cell barcode and UMI tags (the sequences for these are in R2).  Maybe all of the 'cell reads' are in one set of FASTQ files.

Some tools will do alignment *and* quantification.

Kallisto is a pseudo-aligner: it doesn't map to 'reads' to a reference, instead it maps 'k-mers' to a reference transcriptome (a substring, which because it is a window smaller than the gene sequence (read), is not unique to the read, but represents a set of shifted sub-strings).

Bray (2016) says that alignment for k-mers is much faster.  They may also cope better with some sequencing errors.

THe referencing to a transcriptome (rather than a genome, as STAR does), 

"We do kallisto bus on the fastqs of single cells to generate the BUS file and then use BUStools on the generated bus files to get a gene level quantification."

After several steps, the BUStools output will look similar to the CellRanger matrices, namely:

genes.barcodes.txt
genes.genes.txt
genes.mtx


