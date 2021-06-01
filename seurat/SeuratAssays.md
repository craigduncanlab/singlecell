---
title: "Seurat_Assays"
author: "Craig Duncan"
date: "06/04/2021"
output: html_document
---

# Assays

## The nCount and nFeature columns in the main data.frame

The other function that the Seurat object performs automatically is to populate the main data.frame with information about how many counts there are for each given cell, and also how many different genes show at least a count of '1' for that cell.  These are stored in the data.frame using this syntax: nCount_assayname and nFeature_assayname, where the 'assayname' is the currently selected, or default matrix in the current/default Assay object.

When a Seurat object is created, it must contain a default Assay object, which may have one or more of the 'counts', 'data' and 'scale.data' matrices inside it, each of which will be a gene-v-cell count matrix, with genes as rownames.   The Seurat object automatically takes this information for the current Assay object, and does the calculation of total counts detected for each cell, regardless of the genes found.  It also counts the number of features (genes) present in the cell, which is the number of different types of genes, not the toal count.  

For each Assay added to the Seurat object, this cell count information is converted into numeric or integer vectors and added to the main data.frame's columns.   There will be a pair of these columns (nCount_assayname and nFeature_assayname) for each Assay added.

It is worth noting that over time Seurat has moved from a program that was primarily concerned with sequencing of RNA transcriptome data and 'genes' to a more general concern with 'features' (which are wider, but still inclusive of genes).  By Seurat version 4, some of the previous data.frame column references like 'nGenes' have now been renamed to 'nFeature', and 'nUMI' has been renamed to 'nCount'.  In Dave Tang's [Intro](https://davetang.org/muse/2017/08/01/getting-started-seurat/), written in 2017, the feature was still called 'nGene' and the count 'nUMI'.  The data set is still the same, so it's just a matter of getting used to using the new references to follow along his tute.

# Special considerations for Multi-Assay projects

## Seurat assumptions about Assay sizes and count matrices

The insertion of data into a Seurat object relies on the main data.frame containing a list of all the cellnames which is consistent with the information in the Assay objects.   This means users cannot unilaterally change the data.frame or have differences between different Assay colnames or rownames for the same project.   If you want to create a smaller 'subset' of data to work from, you need to create an entirely new Seurat object for that purpose. 

The expectation is that a single 'Seurat' project will update the researcher's work status (and this includes storing classes of cells that are being identified in different ways).  

The Seurat program expects that each data collection will consist of several R Matrices, the second and following of which slightly transformed version of the original data.  To encapsulate these progressive versions of the original count data, Seurat creates the Class of object called an "Assay".  In line with the expected workflows, the Assay object has slots for Matrices that are called 'counts', 'scale.data' and 'data'.

Also, each one of these data collections is usually accompanied by R code that will explain how the project will proceed through several work stages (e.g. raw data in 'counts' is then scaled to create the 'scale.data' matrix, normalised and so on).   There is an expectation that the researcher will, from time to time, change the 'active' version of the matrix in the Assay, and the active Assay that is being used.  This enables the researcher to retain the intermediate data (in Matrices and updating code) to enable reproducible, auditable research.

## Dynamic state, updating functions of Seurat objects

With more than one 'Assay' inside a Seurat object, we would expect that changing the 'active assay' would require the object to also update its meta.data information. {can this be verified?}

In the design of these inter-related Seurat objects, there is conscious attention to the fact that :
1. There may be one or more versions of the raw count data stored in each "Assay"
2. There may be classes of cells that are sets of cell names (defined in an 'Identity' class which stores a range of cells and the identity name).
3. The classes of cells that are created or renamed should be capable of applying to:
(a) one or more of the Assay sets; and 
(b) a particular version of the data (Matrices) within each of the Assays.
4. If there is more than one Assay, the researcher may wish to change the Assay set that has the current focus, by setting the Seurat project object's 'active Assay'.