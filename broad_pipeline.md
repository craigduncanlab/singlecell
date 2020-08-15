# Single-cell analysis – pre-processing for Seurat using Broad Institute pipeline (Drop-Seq/Picard/STAR)

(c) Craig Duncan, April 2020

This picard/Drop-Seq workflow can be broken down into several discrete stages, each with its own sub-steps.   All of these are necessary before being able to use Seurat.

These include:

A.	Obtaining a queryname-sorted BAM

•	Preparing FASTQ files as queryname-sorted BAM

[Picard Tools]

B.	Modifications to annotation files, and creating a genome index using them

•	Modifying GTF files and creating a dictionary (SAM file) (*)
•	Using the modified GTF files to create a refFlat file for Picard Tools.
•	Using the modified GTF files to create the genome index using STAR.

[Drop-seq, Picard Tools and STAR]

(*) If the GTF lacks attribute fields such as gene_name, transcript_name or has a ‘;’ in the gene_name for homologs, it also needs cleaning up ready for the Picard tools.

C.	Preparing aligned data ready for counting

•	Associating cell-barcodes with data in the queryname-sorted BAM to obtain a new, coded BAM.
•	Aligning BAM using the STAR-generated genome index (but losing barcodes information in process)
•	Merging aligned BAM with coded BAM 
•	Generating per-cell counts

[Drop-seq, Picard Tools and STAR]

D.	Output of cell-count ready expression matrices.


There are some published pipelines slightly different to the pipeline assumed in the Dropseq manual.  See for example: https://github.com/mccrowjp/Dropseq (MIT Licence).  Here is another more modified pipeline based on Dropseq: https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-scrna-seq.html#introduction 

Create metadata (including STAR index) for use in later pipeline

Only steps marked (*) below are essential to the Drop-seq/Picard pipeline.  

The rest are considered useful by Broad Institute (e.g. filtering gene_biotypes, making human readable reduced GTF format, and creating the intervals files).

| Package | Function | Input | Output |
| :------------- | :----------: | :-----------: | :-----------: |
| picard | CreateSequenceDictionary(*) | |
| | List of contigs in the aligned FASTQ file (FASTA) and their length. | Reference.FASTA or fasta.gz | Reference.dict (a SAM file)
| | The DropSeq manual refers to fastq file but .fasta or .fasta.gz reference is required (from libraries)	| 	Output SAM file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension |
| dropseq	| FilterGtf |	1. GTF or refFlat, 2. Sequence dictionary (SAM) – optional. | 	GTF (filtered) – without unwanted gene biotypes |
| | A program that will filter any unwanted biotypes from the GTF reference.
| | By including the dictionary (SAM), it also discards any chromosomes that are not common to the annotations file (GTD) and the SAM file.  the G=gene_biotype option defines the gene_biotypes found in the GTF file that should be excluded. | | |
| dropseq | ConvertToRefFlat (*) | Required because the later workflow involves Picard tools, and this requires refFlat file not GTF.  Only annotations found in both the GTF and the SAM file will be retained in the output.	 | 1. GTF or refFlat, 2. Sequence dictionary (SAM)| refFlat file |
| dropseq	| ReduceGTF.  It will only retain features in both the GTF and the SAM file. The default setting is to retain only 3rd-column GTF features that are exon, transcript, or gene. nb: The supplied workflow by Broad institute does not require this to be run before FilterGtf and the creation of the refFlat file.  It is only run before the creation of the intervals files.	| 1. Filtered GTF, 2. Sequence dictionary (SAM) |	Reduced GTF (GTF) |
| dropseq	| CreateIntervalsFiles (*) This is needed for later Dropseq steps… |	1. Reduced GTF, 2. Sequence dictionary (SAM) | Intervals files ( refname.genes.intervals, refname.exons.intervals, refname.rRNA.intervals) |
| STAR	| This uses the basic index mode of the STAR tool to create a genome index.  It uses the original .fasta file referred to above, but it also takes as an input the GTF modified by the filter step (not the reduce GTF step).  The GTF used must be the same one used for creating the refFlat file for use with Picard tools.	| 1. .fasta reference file 2. Filtered GTF	| Genome/index in the STAR output directory. | Includes files:genome,SA,SAindex, genomeParameters.txt |
| bgzip | | |
| samtools | | Zipped fasta? |

A note on Drop-seq programs

The Drop-seq ‘programs’ are, in fact, shell scripts that call on the main Dropseq .jar file with certain program names passed as parameters.  The relevant line in the shell script is:
java -Xmx${xmx} -jar $jar_deploy_dir/dropseq.jar $progname -h

Running Drop-seq programs requires running the relevant shell script that is inside the main Dropseq folder.  So, for example, if you run ./ReduceGTF (without any parameters) inside the main Dropseq folder, it runs the relevant java program.  The output of that program (without specifying any parameters) provides some useful information that explains the motivation of the programmers in this particular case.  E.g.

```
“./ReduceGtf
ERROR: Option 'SEQUENCE_DICTIONARY' is required.

USAGE: ReduceGtf [options]

GTF files are annoyingly complex with a poor definition of what data is in them. So hey, why not write a file parser.
This program reduces the GTF file in to a simplier, easier to parse format, while simultaneously allowing for data to be
filtered.
Version: 2.3.0(34e6572_1555443285)”
```

In effect, this is a stand-alone GTF file parser which has several options (displayed when you type the command above).  The minimum information you need to use it is that you need 2 inputs: the GTF file you want to reduce/filter and the sequence directory (i.e. SAM file).

If this file parser does not work, the resulting GTF file may not be ready for the subsequent workflow, and you will need to take over this role yourself.

```
./ ./ConvertToRefFlat

ANNOTATIONS_FILE=File         The annotations set to use to label the read.  This can be a GTF or a refFlat file. 
                              Required. 

SEQUENCE_DICTIONARY=File
SD=File                       The reference sequence dictionary.  Only chromosomes found in this file AND the
                              annotations file will be retained.  Required.

./FilterGtf

GTF=File                      Input GTF file to be filtered.  Required. 

OUTPUT=File
O=File                        The output filtered GTF file  Required. 

UNDESIRED_GENE_BIOTYPE=String
G=String                      gene_biotype value that flags a GTF record as undesired  Default value: null. This option
                              may be specified 0 or more times. 

SEQUENCE_DICTIONARY=File
SD=File                       The reference sequence dictionary. If specified, GTF records for sequences not in the
                              dictionary will be discarded.  Default value: null.
```

Refs:
1.	https://gatk.broadinstitute.org/hc/en-us/articles/360036712531-CreateSequenceDictionary-Picard-
2.	Source for DropSeq: https://github.com/broadinstitute/Drop-seq/tree/master/src/java/org/broadinstitute/dropseqrna/utils
3.	https://github.com/broadinstitute/Drop-seq/blob/master/src/tests/java/org/broadinstitute/dropseqrna/annotation/FilterGtfTest.java
4.	Issues with GTF formats (particularly when converted from GFF), and need for first column name to match the first column in the reference .fasta file supplied: https://github.com/broadinstitute/Drop-seq/issues/78 

# A note on Picard tools

These are individual java programs contained in the java archive, supplied in the ‘3rdParty/picard’ folder of the Drop-seq package.

To obtain a list of tools and a brief description, just type this from within that folder,e.g.:

```
/Drop-seq_tools-2.3.0/3rdParty/picard$ java -jar picard.jar 
```

To obtain a list of instructions for a specific tool (very detailed, including option parameters, just type the program name as well, without any options or files as input, e.g.:

```
/Drop-seq_tools-2.3.0/3rdParty/picard$ java -jar picard.jar FastqToSam 
```

The information obtained is more detailed than the current Drop-seq manual, and explains more comprehensively the inputs and outputs to each program. 

# A note on the STAR program (this stage: for index creation)

The use of STAR in this pipeline should be created in a *shell file* that sequences the use of the tools, and the required parameters. 

If the STAR index creation process works, you should have these files (or similar), in the dropseq/STAR directory, including genomeParameters.txt as well as the Genome and SA files:

```
Genome  
SAindex
chrName.txt 
chrStart.txt 
exonInfo.tab
genomeParameters.txt
sjdbList.fromGTF.out.tab
transcriptInfo.tab
SA
chrLength.txt 
chrNameLength.txt
exonGeTrInfo.tab
geneInfo.tab
sjdbInfo.txt
sjdbList.out.tab
```

# Pre-processing of the FASTQ data for this pipeline

Preliminary steps in relation to the read data for this pipeline include the need to create a BAM file for the read data to use as an input for STAR that will also be suitable for per-cell analysis.  This will depend on whether you have the raw files (BCL), or a FASTQ file.  

The Drop-seq/Picard workflow requires you have a BAM file and that it is sorted correctly (by first column i.e. queryname).   Steps are needed to convert BCL or FASTQ to SAM or BAM.  

## BCL files transform steps

| | Picard	| Input	| Output |
|:--------|:---------|:-------|:------|
| | ExtractIlluminaBarcodes |	Raw Illumina sequencing data (BCL file)	| BARCODES.Dir for each lane. |
| Demultiplex and sort reads to produce unaligned BAM | IlluminaBasecallsToSam {Not yet FASTQ) | Raw Illumina sequencing data (BCL file) + BARCODES.Dir for each lane + LIBRARY_PARAMS file | Query-name sorted unmapped SAM/BAM (uBAM). |

## FASTQ files transform sequence

Picard provides the ‘SortSam’ tool with the ‘queryname’ option to sort a SAM/BAM primarily according to column 1 (queryname).  

| | Picard	 |	Input	| Output |
|:--------|:---------|:-------|:------|
| Change FASTQ file (paired read) to BAM | FastqtoSam F1=forwardread.fastq F2=reverseread.fastq O=output_unmapped.bam SM=sample_name |  FASTQ (optionally gzipped) | Unmapped SAM/BAM (uBAM) | 
| | SortSam i.e. ```picard SortSam SORT_ORDER=queryname``` | Unaligned, unsorted BAM | BAM file sorted by queryname e.g. a single file called ‘unmapped_bam.bam’ |
| | If you cannot simplify your launch command, you may need: ``` java -jar picard.jar SortSam I=inputfile.bam O=outputfile.bam SORT_ORDER=queryname ``` | | |

# SAM/BAM

The term ‘queryname’ or ‘qname’ is a basic element of the SAM format and is basically the name of each read (https://genome.sph.umich.edu/wiki/SAM).  It is the information in the first column of the tab-delimited SAM format.  Spec: http://samtools.github.io/hts-specs/SAMv1.pdf  (this last spec also has useful information on how to convert from GFF3 to SAM).

A SAM file can be sorted in a different order, and the individual file will have its “SO” tag set to describe the sorting order used.  It can be ‘co-ordinate’, queryname, unknown or unsorted.   The queryname sort will have the result of grouping together duplicate reads with the same name (that happened to be aligned in different place).  Nb there are other fields that contain other useful metadata, including the platform that produced the sequencing data (e.g. Platform/technology used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO).

The SAM spec document contains this important warning regarding the ‘queryname’ sort order: “It is known that widely used software libraries have differing definitions of the queryname sort order, meaning care should be taken when operating on multiple files of varying provenance. Tools may wish to use the sub-sort field to explicitly distinguish between natural and lexicographical ordering. See Section 1.3.1. ”  The natural order is by alphabetic/numerical sequence.

Examples of qname (from the umich.edu site above) include read refs like ```“1:497:R:-272+13M17D24M”.```  The qname need not be unique.

```samtools sort -n``` will be similar (if not same) as the ```sortsam``` in Picard.

# SAM/BAM sort order options

The SAM format provides for different types of sort options:

  “For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field using the reference sequence dictionary (@SQ tag). Alignments within these subgroups are secondarily sorted using the left-most mapping position of the read (POS). Subsequent to this sorting scheme, alignments are listed arbitrarily.

  For queryname-sorted alignments, the tool orders [the] records deterministically by queryname field followed by record strand orientation flag, primary record flag, and secondary alignment flag. This ordering may change in future versions.” 
  
  (https://gatk.broadinstitute.org/hc/en-us/articles/360037594291-SortSam-Picard-)



