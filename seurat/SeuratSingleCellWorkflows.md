---
title: "Seurat_SingleCellWorkflows"
author: "Craig Duncan"
date: "06/04/2021"
output: html_document
---

# Scientific and reproducibility goals

The predominant goal of a Seurat object is for scientific research related to 'single-cell-gene expression-data analysis'.   This means it has to practically store cell data and meta-data, and provide functions to users when they are performing statistical analysis from data input through to visualisation.  This includes:
1. Manage R objects holding gene-and-cell count data called 'Assays', each of which may hold matrices with thousands of genes and thousands of cells. 
2. Allow the user to classify individual cells by category labels (statistical, dimension reduction).
3. To facilitate visual inspection and analysis of 'clustering'. 
4. To manage or use groups of 'markers' for classification of unknown cell types i.e. typical genes that identify cell types, or typical patterns of gene expression that identify cell types.

As it states in this simpler blogs/vigenettes for training:

[UCDavis_Seurat](https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html)

# Single cell sequencing goals

Chen (2016) in BMC Genomics 2016, 17(Suppl 7):508 DOI 10.1186/s12864-016-2897-6 [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5001205/pdf/12864_2016_Article_2897.pdf) said:
"Unlike the objectives of conventional RNA-seq where differential expression analysis and the detection of differentially expressed genes (DEGs) are integral components, the most important goal of scRNA-seq is to identify variably expressed genes (VEGs) across a population of cells to account for the discrete nature of single-cell gene expression and uniqueness of sequencing library preparation protocol for single-cell sequencing.

# Wet lab and experimental aspects of RNA seq

The scientific goal is touched on in Macosko, at page 1212:
"Ascertaining transcriptional variation across individual cells is a
valuable way of learning about complex tissues and functional
responses"

The 2015 Macosko study into gene expression used Mouse and Human cells (creating barcoded libraries) was an early study illustrating Seurat-based analysis in the context of drop-seq methods for sequencing.  *This paper, more than any later applications  of this approach in the plant field (like Shulse etc, pointed to by those who have no better idea), provides insights into the significance and purpose of the statistical analysis which is inherent in the Seurat program's functions*  The article uses the acronym 'STAMPS' to refer to all single-cell transcriptomes (identified with barcdoes) attached to microparticles (the raw drop-seq library).  They are distinguishable from beads never
exposed to cell lysate.  

Macosko also identifies some potential features in single cell drop-seq lab processes:
1. Potential for cells to stick together (doublets)
2. Ambient RNA levels high due to cells damaged at start of process.
3. PCR and sequencing errors that inflate the apparent unique UMIs (so collapse these - see page 29/51 in Supplementary Extended information with the Cell paper)

In early bulk sequencing of RNA [ref](doi: 10.1093/bioinformatics/btw174) and in single cell (see Macosko et al. 2015), the majority of STAMPS were exposed to ambient RNA (i.e. not even in the experiment's cells), so the very low RNA 'read' counts in many cells (e.g.>90% with reads <500) are eliminated (cut-off point 'knee'), to leave the data with only the useful (nucleic) RNA data.  This is usually done by arranging the cumulative number of counts per cell from lowest to highest, and where there is an abrupt 'knee' and levelling off, you can eliminate all of the lower read STAMPS.

Is the pbmc data pre or post this?  We can have a look at the data to see.  In Dave Tang's 2017 [intro](https://davetang.org/muse/2017/08/01/getting-started-seurat/), he noted that there were 69,000 'reads' per cell, which seems quite high, but apparently nowhere near as high as bulk RNA sequencing.

The tutorial goes through and eliminates cells based on those with low numbers of detected genes (which will probably correlate to low absolute 'reads' as well, but not necessarily).  This is also influenced by the ability of the alignment/annotation process to identify a known gene.  That is, gene expression analysis is always limited by the known genes that are used for the statistical comparisons.

## OSCA (Bioconductor) comparison

It's interesting that the Macosko data is referred to in the OSCA (Bioconductor) quick start, and is even included in the 'scRNAseq' (R) package:
[OSCA5.5](http://bioconductor.org/books/release/OSCA/overview.html#quick-start)

There they follow a generally-recognisable gene expression data analysis workflow.  This includes filtering, normalisation, performing some PCA (initial 'low rank representation') and using t-SNE or nearest-neighbour clustering, before visualising [expression] clusters in a UMAP plot.  Seurat follows most of these steps, but the data is held in a different R object (the Bioconductor workflow there holds the data in a 'singleCellExperiment' object). 

The Quickstart page has useful advice, but the specifics of how to prepare count matrices require experience, and practice in problem-solving using actual data.

# Getting to counts and beyond

A useful workflow diagram in rCASC paper (Figure 1) in (2019) [GigaScience](
https://academic.oup.com/gigascience/article/8/9/giz105/5565135)

Genomics 10X and DropSeq UMI are represented as being different workflows to get you to a counts matrix.  (Presumably DropSeq one is without CellRanger).  Genomicx 10X 'freely' provided a pbmc dataset (10X) which is used by Seurat for tutorials, since it has already been processed up to the point of CellRanger output, where it can be read in (these are not zipped files in the pbmc demo).  The dataset is 'filtered' which as XGenomics explains, "Contains only detected cellular barcodes" (i.e. not based off a whitelist).  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

TO DO: Read http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression

Also, interestingly, I've found that the CellRanger software has opensource for compilation (3 years old?).  It is found here, along with a link to downloadable binaries (subject to licensing, non-commercial use etc): https://github.com/10XGenomics/cellranger

# Scientific knowledge about single cell analysis needed to use Seurat

In order to arrive at the use of Seurat with the minimum information needed to understand this area of scientific research, some familiarity with these topics will be needed:

1. The relationship between cell type evolution (specialisation over the life cycle), the theory that gene expression (the number of genes expressed, and therefore counted in the sequencing process) is different in different cell types.
2. The ability to capture gene expression levels by reading the transcriptome (RNA seq, single cells).
3. The production of a counts matrix (genes expressed versus cells, and counts).  Using STAR, Hi-Seq, FeaturesCount.  (It is helpful to know that in this context, genes are often called 'features', but in other contexts, features can refer to other things).  It is the fact that the sequencing process fragments cells into molecules, which will each contain different fragments that may contain different gene sequences, that requires a deliberate 'accounting' process to match particular genes to particular fragments, and then to collate this information for all the different cells found.   There is a need for a reference index to what gene labels can be applied to different gene sequences (transcription library file), and a procedure for recognising the barcodes associated with reading of fragments.
4.  The idea that expression of certain genes can be used as 'markers' of certain cell types, in that they tend to vary by cell type (providing points of differentiation).  This is associated with the idea that classification can be done using the evidence of the presence/absence of genes (the available data). so gene expression can be interpreted as a 'signature' (marker) of cell type.  When this is well established these are 'canonical marker genes'.
5. A whole host of statistical techniques might be brought to bear on the question of what 'differences' in gene expression can be correlated with known or unknown gene expression types. This includes dimension reduction, tSNE, PCA.
6. Investigation techniques or ways of determining relevant genes are assisted by data visualisation techniques (which includes clustering, mapping onto a 'UMAP').   Depending on the state of the work, the presence of common genes in the transcriptome can be used to associate expression levels with cell types, and visualise them as 'clusters'.  Once these are identified and correlated, there can be labels attached to the clusters, either to identify the genes used for classification, or the presumed cell types associated with these clusters.
7. The focus of the study can emphasise one or more of these elements (e.g. models for classifying cell type by gene expression, machine learning algorithms and training, or just identifying unknown cell types, or identifying new gene expression markers for cells).   To some extent, Seurat's tools can make it easier to proceed to a simple 'identification' of cell types present by choosing conventional statistical, clustering and marker sets in order to get a basic clustering and identification.  The user may still need to manually annotate some of the cells identified visually in a cluster, by creating labels and sets.

The skills required for all of the above include scientific knowledge, scientific research techniques, statistical knowledge, computing knowledge and, in some cases, programming ability.

## Seurat's role

A tool like Seurat is customised for performing the statistical transformations and analyses easily, from conversion of data through to visualisation and retrospective interpretation. The authors of Seurat have made choices based on current statistical techniques, and, in particular, choosing those like t-SNE that seem to provide more useful data visualisations.

The net result is that Seurat provides a set of functions that summarise a number of complicated statistically-based data analysis techniques or goals in just a few lines of code.  They are really the tip of the iceberg in terms of what is being performed mathematically, but that is the point.   A user of the system needs to be able to summarise the complexity of the statistical steps and data analysis, in order to wisely leverage the power of the work invested in the development of Seurat.  

The purpose of this review is to offer some high-level comments on the Seurat package, derived from my inspection and experimentation with the package and R language. The 'Seurat object' is the major custom R 'object' provided to users as part of the package.  

This manual is ultimately intended to help with single cell analysis at the point in the workflow at which Seurat becomes applicable.  I have also made some introductory comments about the earlier stages in the sequencing process, because it helps explain the inputs to Seurat objects (gene-vs-cell counts), including how they are prepared and the extent to which quality controls should have been implemented prior to that stage of analysis.  

As explained further below, this manual is not designed to merely teach you how to use Seurat, but how to practice as a bioinformatician, and to carry out single cell analysis in the context of being a scientist.

You can directly compare the workflow for the pbmc data discussed in the Seurat tutorials with the same data in the OSCA/[Bioconductor workflow](http://bioconductor.org/books/release/OSCA/filtered-human-pbmcs-10x-genomics.html), that uses libraries such as scater, scran and a specific 10X data loading package called TENxPBCMData.

# Pre-requisite knowledge/knowledge broadening

To master working on single cell genomics analysis with Seurat, users or teams will need to have some or all of these skills (a program for grad students working in bioinformatics):

- biological knowledge of how sequencing is conducted, and how this affects data, errors in sequencing data collection 
- detailed knowldge of formats for sequencing data (digital)
- detailed knowldge of gene libraries for common genomes (plant, animal)
- detailed knowledge of the preparation of gene-v-cell-count matrices, and data formats for this
- detailed knowledge of statistics, probability distributions, dimension reduction (PCA, t-SNE) and clustering (for data science in general)
- knowledge of R language and base data structures (for scripting and also understanding Seurat package itself)
- understanding of how Seurat object is custom data object
- understanding of functions provided with Seurat package work with normal scientific research/data analysis workflow
- knowledge of marker genes that will be useful for retrospective identification of cell types

Bioinformatics presents a specific (but not unique) problem in that it requires sophisticated solutions to problems of working with large data sets and memory limitations or speed on computation hardware.  It is therefore advisable that any scientist who happens to be using Seurat should also learn the high-level features of the process and whether they can be implemented in a different hardware or software ecosystem, possibly (but not necessarily) with some trade-offs in terms of speed and precision. For example, the workflow in Satija lab's Seurat tutorials for single cell analysis from reading in .mtx data through to clustering can be substantially repeated in a python language environment using pandas and scanpy packages, as illustrated by the single cell sequencing Clustering Tutorial for pbcm10k data, in this [webnotebook](https://crc.pitt.edu/sites/default/files/scrnaseq/pbmc10k.html).

# When to learn/use Seurat?

In my view, because Seurat is an all-inclusive package designed to it is best to use and familiarise yourself with Seurat *after* you are already familiar with the biological and statistical procedures it is designed to present to the user in summary form.   These processes should look familiar enough to you that you understand the difficulty that Seurat is attempting to overcome, and you understand why it is a labour-saving choice.   You should understand the alternatives to it (both in terms of the more laborious, mathematical techniques and also the software prepared in competition with it). (If you arrive at Seurat first, you have to go about this backwards: understanding the context for it, then returning to it with fresh eyes).

By all means, use Seurat to speed up the process of exploring and comparing data-sets.  However, you can't usefully do this (or do it quickly and accurately) until you have the underlying knowledge.  You do not 'learn single cell analysis' or 'statistics' or even 'machine learning' merely by following Seurat tutorials. In fact, Seurat is really just a high-level shortcut for those who already understand, broadly, what single cell analysis is, and why you would want (need) tools to provide statistical analysis and visualisation.  A tool like Seurat lowers the bar for how much work you do, and might even offer some possibilities of non-specialists getting the same outputs, but doesn't really lower the bar for the general knowledge you need for single cell analysis.  Statistical analysis can provide a way of selecting, reasonably, subsets of the data set, and Seurat is merely a tool to help take out some of the labour.  

Do not expect to use Seurat to quickly enable you to do an analysis and write a paper (on your first use or even your 10th use of it).  The more experience you have in exploring different data sets with similar processes (and testing, observing natural variation), the more insights you can bring to any particular data sets you are investigating.   Alternatively, if you are merely using Seurat to help someone 'get things done', then you will need to be supervised by a person with the pre-requisite knowledge described above.
