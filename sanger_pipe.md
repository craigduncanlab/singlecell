# Sanger Institute Pipeline 

(c) Craig Duncan 2020

The Sanger institute, associated with the development of the FASTQ format, states: 

  “We have publicly available perl scripts capable of demultiplexing any scRNASeq data with a single cell-barcode with or without UMIs for plate-based protocols.” 
  (see Sanger Ref 1 and 3) 

  The relevant readme/perl files for this purpose were last updated 2016-2018 and include perl files to gather summary statistics, including counting reads, and obtaining the number of unique ells and the number of reads for each one.  The files include: 

```
00_Steps (summary of the shell scripts to be run, and why) 
0_Determine_Barcodes.pl 
0_Check_Barcodes.pl 
0_Gather_Summary_Statistics.pl 
```

Sanger refs: 

1.  https://scrnaseq-course.cog.sanger.ac.uk/website/processing-raw-scrna-seq-data.html 

2.  https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_pre-QC.htm 

3.  https://github.com/tallulandrews/scRNASeqPipeline  

4.  https://scrnaseq-course.cog.sanger.ac.uk/website/construction-of-expression-matrix.html 
 
