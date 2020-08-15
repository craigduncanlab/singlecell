# The demultiplexing problem for raw scData 

The term ‘demultiplexing’ is often used to describe the general process of separating mixed data streams into separate output files.  Some sequence service providers describe demultiplexing as the process of providing lane-based FASTQ files from a single sequencing 'base call’ file (BCL).  See Demultiplex Ref 1.  However, in relation to single-cell analysis, it refers to the identification and possible separation of cell-specific information from within each of the FASTQ files. 

There is a necessary and important data-separation and identification step involved in the sequencing process.  For one or more samples (and there may be many cells each constituting a sample), there are usually two FASTQ data files (reads) prepared, in which one has the barcodes for a cell, and a read reference, and the other has a read reference and the biological data.   It is the biological data that is aligned to the reference genomes, so there is a risk that this is done but there is no retention of information about which cell was involved. 

Since cell analysis depends on the prior step of using barcodes to mark cells in the sequencing process (ligand bonded), to identify each individual cell uniquely for each read.  For massive sequencing projects, a single (or few) FASTQ file includes read and barcoding information for reads of many thousands of individual cell samples, but they are not usually grouped by cells.  The goal of later identifying the relationship between individual cell samples and reads is called ‘demultiplexing’.  This demultiplexing step (probably better called the cell data grouping) is essential to preparing the data for cell-based analysis (where observation data is on a per-cell basis, rather than in bulk).   

How the cell information is identified with specific biological information in reads depends on experimental choices and available sequencing information. Some sequencing providers provide their own software to achieve these goals (for example 10x Genomics has CellRanger), but not all researchers wish to use this, given licensing restrictions.   

If you are not willing or able to use CellRanger, then there is a preliminary question about how you can obtain and retain cell-specific read data (and ‘de-multiplex it’) at an early stage, and also preserve it after the alignment step.    This threshold data problem has resulted in several different workflows and ‘work-arounds’. Demultiplexing often involves custom scripts to suit the protocol and the workflow. 

The use of ‘whitelists’ (legitimate barcodes for the samples) is another issue.  If sequencing providers have these, they may not freely publish them.   A choice needs to be made between using those programs that require a whitelist, and those that do not.  A further complication is whether to generate a whitelist, just to be able to use a program that requires it. 

How to deal with the cell-specific biological information (obtained with reference to barcodes) is another matter of experimental choice.  One choice is to filter the original master FASTQ file and create distinct output files, one FASTQ file per cell sample.  Separate files may allow for some parallel processing advantages.   It is at least possible to split the data into separate files, using a tool like Sabre (see Demultiplex Ref 1), which also requires a separate whitelist (reference) file for the barcodes (see Demultiplex Ref 2).  Sabre seems fairly simple in that it does not seem to allow setting any minimum read count depth before outputting individual .fq files.   

Since the original data is often supplied with barcodes in one FASTQ file, and biological data in the other, another approach is to re-organise the barcode information in the few FASTQ files that are given, so that it is stored with the biological sequencing data, either before alignment or somewhere else in the pipeline (see Sanger Ref 4 below at para 4.6.2).   The UMI tool assists with this but does not demultiplex barcodes to get one fastq per-cell. Instead reads are tagged with the cell barcode but kept together in one file.   A common approach is to transfer the cell identity data (barcodes) in the barcode (Read1) file to the header section of the biological data file (Read2).  See Sanger Institute and Broad Institute pipelines below. 

Whichever choice is made will influence the subsequent workflow to achieve the data analysis goals. 

Demultiplex Refs:

1.	https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/mkfastq

2.	https://astrobiomike.github.io/amplicon/demultiplexing 

3.	https://github.com/najoshi/sabre 
 
# Whitelists 

The use of whitelists is another form of quality control and verification.  Barcodes in FASTQ files do not necessarily correspond to an actual read, and some are based on other material in the cell to the desired target.   Many of these other possible experimental markers will show up as small reads.  If possible, when attempting to verify that barcodes found in a FASTQ file are legitimate.  a master list of legitimate barcodes is used, in combination with a barcode search in the FASTQ file, to confirm that the barcodes found are not invalid.  This ‘whitelist’ is sometimes a proprietary list, so it may not be available if not supplied for some reason.   

The STARsolo module in STAR, which attempts to assist with extracting gene expression data, still requires a whitelist, apparently with the same content as CellRanger’s whitelist requirement. 

If you don’t have a whitelist, can you make one that is reliable enough to be used?  Probability-based assessments of which barcodes that appear are likely to be valid are an alternative approach (even used in the CellRanger software pipeline itself).  Some tools like UMI tools and buslist seem to apply this approach.  UMI tools is a python program that will accept a FASTQ file and output a predictive whitelist.  Drop-Seq workflow, which does not require an explicit whitelist, still has elements that require the user to form decisions about threshold tests for what is a legitimate barcode.  These will involve similar considerations about the frequency of barcode detection which are similar to what the independent whitelist-creation tools use. 

Whitelist Refs: 

1.  https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist- 

2.  https://umi-tools.readthedocs.io/en/latest/reference/whitelist.html 

3.  https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html 
 
