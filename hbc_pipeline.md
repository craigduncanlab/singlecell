# Harvard Chan Bioinformatics Core (HBC) pipeline 

(c) Craig Duncan 2020

The approach taken to demultiplexing and handling cell barcodes seems to be to take these from the Read 1 file and put them into the Read 2 file header (see Harvard Ref 1).  The python program that can take all non-molecular data (both cell barcodes and UMI) is called ‘umis.py’.  It is found at the Harvard Ref 1. It is run with the command umis fastqtransform. 

One tutorial suggests limiting the analysis to cells by filtering those with <1000 reads (Harver Ref 2) 

Harvard Refs 

1.  https://github.com/vals/umis and https://github.com/vals/umis/blob/master/umis/umis.py  

2.  https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_pre-QC.html 
