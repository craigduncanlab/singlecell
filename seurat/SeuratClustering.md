---
title: "Seurat_Clustering"
author: "Craig Duncan"
date: "06/04/2021"
output: html_document
---

# Clustering

## Clustering is a general scientific method and concern

The central importance of 'identities' to a Seurat object reflects the scientific concern with annotating, measuring or classifying cells and cell characteristics by different schemes (cell type, quantitative markers, quality cut-offs).   

The usual situation is that a Seurat project will end up with many new additions to its main data.frame, to reflect new classifications useful to analysis and visualisation.  This classification and attribution may, in fact, be one of the most important outcomes of the work.

See this:
[ClusterLink](https://www.r-bloggers.com/2020/05/how-to-determine-the-number-of-clusters-for-k-means-in-r/)
Questions: what is data set?  How do we identify inherent dependent variables?

Notice how clusering, elbow-plots (which are methods identified in bioinformatics, even single cell genomics) are not unique to that area but are:

1. General statistical terms/methods
2. Useful in biological in general
3. Still useful in relation to genetic analysis
4. Still useful in relation to the data that is captured in a genetic sequencing analysis.
5. The subject of software tools that are designed for any one of the above 4 areas.

You may find that some custom tools in bioinformatics are, in fact, wrapping pre-existing tools into a workflow or software object, but hiding it for convenience.  

## The use of cell identity schemes and modification of the main data.frame 

Why is there so much attention in the Seurat object's structure to the notion of 'ident', and several different concepts?  

The Seurat object seems to have evolved so that users can specify 'idents=' parameters in some general Seurat functions, in order to choose subsets of the cells.  This in turn requires users to apply a particular scheme of classification to all cells in the project.

To achieve this goal, each new Seurat object has:
1. Automatic creation of its main, default data.frame with all cellnames as rownames
2. Automatic population of column in that data.frame with a classification scheme - 'orig.ident'.
3. Automatic definition of that classification scheme to reflect the project name i.e. 1 level.
4. Automatic setting of active.ident to 'orig.ident'
5. Housekeeping functions for these identities. (DefaultIdent(), Ident() etc)
5  General Seurat functions that take the (idents=c(...)) parameter, to set the *levels* to filter on, within the current identities scheme.

The housekeeping functions are all based around the idea that any kind of cell-category list or scheme must be:
1. Based on (or equivalent to) a column in the main data.frame, which has cellnames as the rownames. 
2. The classification scheme will apply to all cells (even if explicit 'levels' are only given for some cells and some cells have no 'level' attached to them for the active.ident?) 

Though the Seurat object does not require users to choose which data.frame to use for meta-data, it does anticipate some of the columns in that data-frame being selectively used for cell 'identities', or 'clustering groups' or just 'categories for filtering'.   Users will be able to add new data.frame columns, and to choose the active column that will be used for filtering in some of the other Seurat functions.  

To use the identity scheme with Seurat functions effectively, your R code must work your way through the housekeeping functions first before calling the relevant function that has an 'idents=' parameter.

The basic procedure to use idents is like this:
1. Prepare some data column for the main data.frame (initially it's orig.ident which only has 1 level).  This will have some set of 'levels', based on the work to day.
2. Set the active data.frame to the one you want with 'DefaultIdents()' function.
3. At this stage, the identity 'factor' set as default will be the value that R reports is in the object@active.ident slot.  It will have the same length as number of rows with cellnames.
4. Use your function, and choose from the available levels in order to filter the cells using those levels as a subset of the identities for the cells.  i.e. function(....idents=c('PCgroup", "mtGroup") ).

The advantage of adding columns to the data.frame is that the cellnames are already contained in the rows (and these may number thousands of cells).  By manipulation of existing identity columns, users can arrive at new categories for filtering based on the existing data, rather than start from scratch.

Users add data.frames in the way they would normally do in R.  Choosing the active column for the Seurat object's identity functions is done by way of the 'Idents()' function, which is used to specify which of the data.frame columns is used to specify cell clusters.

One of the general assumptions that Seurat makes about what kind of data is in the column being used for 'idents' or 'identities' is that it will be one of R's 'factors'.  This means that there will always be 'levels' implicitly associated with that data.  Users can modify factor levels in base R, if required.

If a data.frame column is set to be the 'identity' column with the 'Idents()' function, then it is internally stored as a 'factor' data type (somewhere), even if the source data.frame column is not a data type.  The Seurat package does not force the original data.frame column to be factor, just the information that is set as the 'Identity'


