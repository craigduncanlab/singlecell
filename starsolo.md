# STAR AND Starsolo 

2019 - Now STAR has been updated, to do some ‘Cell-Ranger’ tasks (but very much based on CellRanger, and the nature of the data CellRanger works with or will output) 

STAR release 2.7.0c (10 Feb 2019) included starSOLO for the first time: https://github.com/alexdobin/STAR/blob/master/RELEASEnotes.md 

starSolo - based on CellRanger counting policies: https://github.com/alexdobin/STAR/issues/737 

There are still some practical issues with using STARsolo, because it doesn’t seem to work totally independently of CellRanger and the information you need to use with CellRanger (e.g. barcode lists).: 

1.  STARsolo asks you to check for number of base pairs in R1 files (e.g. 28bp).  This tells you it is Chemistry option : Cell Ranger v3 
2. You need the CellRanger whitelist files for barcodes (3M-february-2018.txt.gz). i.e. the list of all possible bar codes included in the assay kits. 

CellRanger v 3 has a Read1 with 16bp of CB and 10bp of UMI at the end. (Galaxy Ref 2) 
 
The whitelist comes with Cell Ranger (3.0 or higher).  Where is this available to non CellRanger clients?  *If you don't have this, can you substitute?* 

For discussion see Galaxy refs 1 and 2. 

A suggested workflow (and example at Galaxy Ref 2) is: 
1. Extract Headers and Sequence from a FASTQ Barcode file 
2. Extract 10X barcodes from the sequence file 
3. Rank barcodes by frequency 
4. Subselect a wide range of good to low quality barcodes 
5. Map the barcodes back to the sequences and then back to the header they came from 
6. Filter original FASTQ with wanted headers 
7. Result 

Galaxy Refs 

1.	https://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/scrna-preprocessing-tenx/tutorial.html 
2.	https://zenodo.org/record/3457880/files/subsetting_data.txt 
