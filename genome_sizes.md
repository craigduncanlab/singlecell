# Genome sizes - some draft notes

(c) Craig Duncan 2020

# An illustration

Human genome: 3,257,319,537 bp

bp = 

 - A-T(DNA) or A-U(RNA); or
 - G-C (DNA and RNA)

Gbp-->
Wheat, barley, rye, pea?
| Species | Gbp | Name | Ensembl |
|:-----|:-----|:-----|:-----|
| H.Sapiens | 3.26 | Person | |
| H. vulgare |	4.79	| barley | |
| S. cereale |	6.67	| rye | |
| T. dicoccoides |	10.50 |	Emmer wheat | [link](http://plants.ensembl.org/Triticum_dicoccoides/Info/Index) |
| T. aestivum  |	14.50 | Common wheat | [link](http://plants.ensembl.org/Triticum_aestivum/Info/Index) |
| P. polyphylla (PPY) | 82.55 | Paris | | 

nb useful index of common names and Latin terms [https://www.ncbi.nlm.nih.gov/books/NBK208347/](https://www.ncbi.nlm.nih.gov/books/NBK208347/) and [here](https://www.inspection.gc.ca/plant-health/grains-and-field-crops/list-of-grains-and-field-crops/eng/1323244558875/1323244642996)

# Relative sizes

## DNA

DNA can be much longer (terminated for each chromosome e.g. 5 centrimetres long for each, humans have 46 of them, roughly paired into 23 'sets' of chromosome pairs).

It has been said that human DNA would equal a length of 2.3 metres long in every cell when unravelled.   An example, is [here](https://hypertextbook.com/facts/1998/StevenChen.shtml), where it is said that the average length of the 46 chromosomes in the human body is 5cm (46 x 5cm = 230cm = 2.3m).  In that calculation, it states that there are closer to 6 Gbp per cell (i.e. twice the human genome nt number).

Does this make sense?

2.3m = 2.3 x 10^9 nanometres (nm)

If for human DNA, there are 3.26 Gbp (3.26 x 10^9 nt) in a length of 2.3 x 10^9 nanometres, then this equals only about 1.4 nt per nanometre.

Conventional wisdom [e.g.](http://www.dietzellab.de/goodies/numbers.html) is that the base pairs are 0.34 nm apart (so closer to 3 nt per nanometre)

So why the difference?   The reference to base pairs in the 'human genome' is to the amount of genetic material needed to define the human species in 'haploid' cells (like germ cells: eggs and sperm).  This is the minimum amount of sequential 'code' that describes the genetic material, but that's not actually how DNA is packed into most of the other, 'somatic' cells, which contain two copies of helical DNA, inherited from each of the parents, and this explains our chromosome 'pairs'.  So in somatic cells, there is 2 x 3.26Gb of genetic material and this is the total 'amount' of DNA that is commonly referred to when the length comparisons for ~2.3m are made.

nb: A more recent (2019) academic discussion of DNA metrics is contained here:

Piovesan, A., Pelleri, M.C., Antonaros, F. et al. On the length, weight and GC content of the human genome. BMC Res Notes 12, 106 (2019). https://doi.org/10.1186/s13104-019-4137-z   [link](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-019-4137-z)

## RNA

RNA might only be a few thousand base pairs long (also, it has a ribose (sugar) backbone and is less stable than the 'de-oxyribose' DNA, especially in alkaline conditions; RNA gets attacked by enzymes).

Ribosomal RNA typically represents the largest portion of the genetic code representing an RNA sequence.  In the organisms's DNA sequence itself, the regions coding for the rRNA should be found in specific chromosome(s) and locations.  The codes may be repeated, and there may be several steps between initial translation of coded sequences (sometimes called 'precursor' regions) and the production of smaller units needed for biological activity.

In transcribed RNA, especially those rRNA involved in DNA transcription and protein building, the types are identifiable and predictable in size. They are often found tightly located within 'units' in the ribosomes. 

| Type | RNA Sequence | rRNA type | Size (~) |
|:----------|:----------|:--------|:--------|
| human rRNA | full | rRNA | 7216 nt |
| eukaryotic | 5S | rRNA | 120 nt |
| prokaryotic | 16S | rRNA | 1542 nt |
| eukaryotic | 18S human | rRNA | 1869 nt |

Refs:

[(2018)Epigenetic expression of 5S in A. Thaliana](https://pubmed.ncbi.nlm.nih.gov/29518237/) 
"A typical Arabidopsis 5S rRNA gene is 500 bp long, comprising a 120 bp transcribed sequence and a 380 bp spacer region.":
[PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5887818/)

[5S RNA database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC99124/)

[16S examples](https://www.ncbi.nlm.nih.gov/nuccore?term=33175%5BBioProject%5D+OR+33317%5BBioProject%5D)

[18S5 link](https://www.ncbi.nlm.nih.gov/nuccore/NR_003286.2)

[28S link](https://www.ncbi.nlm.nih.gov/gene?Cmd=DetailsSearch&Term=100008589)

Two units of measurement:

 - nt = nucleotides, which is essentially a length, and 
  - S = Sverdberg units where this refers to settling time in a centrifuge (so S units cannot be added together meaningfully)

mRNA - in eukaryotes, messenger RNA copies (transcribes) portions of DNA (i.e. it is assembled at DNA sites) and transported to the ribosome 'factories'.  Also, the mRNA must also have non-coding regions removed (spliced) in order to become 'mature'.  Unlike the relatively predictable size of rRNA for an organism, the size and length of mRNA varies depending on the genetic material that is copied into the mRNA coding sequence.

tRNA - brings in amino acids so the ribosome factories can make proteins (translation).  Think ''amino acid transfer', 'transport'

rRNA - the RNA that is part of the structure of the ribosome factories: initiates the protein building functions. i.e. you know this is localised to the mRNA.  Due to the persistency and consistency of rRNA across organisms, it is often used for evolutionary studies. 

Refs:

https://www.technologynetworks.com/genomics/lists/what-are-the-key-differences-between-dna-and-rna-296719

## Coronavirus

Compare above ribosomal RNA sizes to:

| Type | RNA Sequence | rRNA type | Size (~) |
|:----------|:----------|:--------|:--------|
| coronavirus | full | RNA | 29903 nt |

SARS-Cov-2 genome size: 29,903 nt 

[posting](https://virological.org/t/novel-2019-coronavirus-genome/319)

[library link](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)

Virus stability issues:
2020, Wakida et al, Stability of RNA sequences derived from the coronavirus genome in
human cells, Biochemical and Biophysical Research Communications 527 (2020) 993e999
https://doi.org/10.1016/j.bbrc.2020.05.008

# Genome size definitions

"Genome size is the total amount of DNA contained within one copy of a single genome. It is typically measured in terms of mass in picograms (trillionths (10−12) of a gram, abbreviated pg) or less frequently in Daltons or as the total number of nucleotide base pairs typically in megabases (millions of base pairs, abbreviated Mb or Mbp). One picogram equals 978 megabases.[1] In diploid organisms, genome size is used interchangeably with the term C-value. An organism's complexity is not directly proportional to its genome size; some single cell organisms have much more DNA than humans (see Junk DNA and C-value enigma). "

Genome size discussion : [P. Polyphylla](https://www.biorxiv.org/content/10.1101/2020.06.01.126920v1)

"Today, the term is used in 2 distinct ways to indicate either the total number of genes or the whole amount of nuclear DNA"
see (2012) Cytogenet Genome Res 2012;137:97–112
DOI: 10.1159/000338820

# Chromosomal explanations for genome size

Assembly Issues

 - Sequencing - entire genome is sequenced many times over
 - Assembly : attempt to assemble connected (contiguous) bases (for all chromosomes)
 - SNPs : For diploid organism, this is differences in single bases between the two chromosomes contributed by the parents (in gametes). *single nucleotide polymorphisms*. 
 - phased assembly: ability to identify gamete sources for all SNPs. 

Examples of differing chromosome sets.  These categories tell us how many replicated sequencies of DNA material we are likely to find.  Not that they are identical in each chromosome set, but that the sequence in each chromosome is comparable and potentially involves the ability to swap out genetic material at the same place in the sequence.

Diploid organisms:
 - H. Sapiens. (i.e. 23 'paired' chromosomes = 6Gbp in total, but say only 3.26Gbp needed to define the human 'genome')
 - Drosopila

Haploids:
 - bee

Polyploids:
 - Wheat
 - Paris Polyphilla (tetraploid; some are diploid)

# Genome size benchmarking

In any case genome size (GS) varies enormously, without apparent 'patterns' ... 

## Mammalian

"It has been found from data-mining GS databases
that GS is a useful cyto-taxonomical instrument at the level
of orders/superorders, providing genomic signatures char-
acterizing Monotremata, Marsupialia, Afrotheria, Xenarthra Laurasiatheria, and Euarchontoglires. A hypothetical ancestral mammalian-like GS of 2.9–3.7 pg has been suggested."


"For more than 50 years, we have been collecting data
on genome size (GS) and have discovered that the ge-
nomes of eukaryotic organisms vary over 200,000-fold in
size."

"When speaking of GS, the histori-
cally based consensus on the meaning of GS should be
clear: the total amount of a mature germ cell’s DNA. The
‘C’ terminology, introduced by Hewson Swift in 1950 in
an unfortunately ambiguous manner (since he used C to
refer to class, category, content, or constant in his first
paper on animal nuclei written for his PhD dissertation),
has been resolved thanks to a personal communication
by H. Swift himself to M.D. Bennett [Greilhuber et al.,
2005]: the term ‘C-value’ was intended to mean ‘constant’
(i.e. the amount of DNA that was characteristic of a par-
ticular genotype) and it was defined as ‘DNA content of
the unreplicated haploid chromosome complement’
[Bennett and Leitch, 2005]."

"In fact, the C-value
enigma conceptually better encapsulates all the relevant
questions related to GS variations. These questions in-
clude: what is the type of noncoding DNA that gives dif-
ferent genomes different sizes? What is the origin and the
intra- and intergenomic transfer of this noncoding DNA?
What, if any, are the functions (at any hierarchical orga-
nization level) that these DNA sequences provide?"

Ref: 

Redi C, A, Capanna E: Genome Size Evolution: Sizing Mammalian Genomes. Cytogenet Genome Res 2012;137:97-112. doi: 10.1159/000338820

[link](https://www.karger.com/Article/Abstract/338820)

Summary :

A hypothetical ancestral mammalian-like GS of 2.9-3.7 pg (picograms)
Monotremata (∼2.97 pg) 
--->Afrotheria (∼5.5 pg) and Xenarthra (∼4.5 pg)
Marsupialia (∼4.07 pg), 
Euarchontoglires (∼3.4 pg)
Laurasiatheria (∼2.8 pg
(546 mammalian species sized from 5,488 living species) 

##  Plants

[P.Polyphilla](https://www.biorxiv.org/content/10.1101/2020.06.01.126920v1)

### Wheat, Barley and Rye

"Any project seeking to deliver a plant or animal reference genome sequence must address the question as to the completeness of the assembly. "

"Seventeen years have passed since the joint announcement of the human genome sequence [14,15]. This period has seen a number of attempts to complete the assembly, applying a variety of technologies [16,17]. All of these have reported a smaller genome size than what has, as of the end of 2017, been suggested in GRCh38.p12, the most recently released Genome Reference Consortium version, which comprises 3,257,319,537 bp. Assuming the Doležel et al. [18] conversion of 1 pg = 0.978 Gbp, 3.5 pg 1C DNA is equivalent to 3,423,000,000 bases. Thus, the 7 pg value represents an ~5.1% over-estimate of the GRCh38.p12 assembly prediction. This difference lies at the lower end of the error range predicted by Doležel and Greilhuber [13]. Given that the human reference genome is still incomplete, the expectation is that the gap between the 7 pg figure and the “real” human genome size will continue to diminish. Nevertheless, a 5% error is not dissimilar to the variation observed between estimates of nuclear DNA amounts of a given species produced by different laboratories [19,20]. Thus, the recommendation remains that the 7 pg figure continue to be used as the reference for measuring 2C values of both animal and plant genomes."

Dolezel, Jaroslav & Janačížková, Janačížková & Simkova, Hana & Bartoš, Jan. (2018). One Major Challenge of Sequencing Large Plant Genomes Is to Know How Big They Really Are. International Journal of Molecular Sciences. 19. 19113554. 10.3390/ijms19113554. 

[link](https://www.researchgate.net/publication/328875915_One_Major_Challenge_of_Sequencing_Large_Plant_Genomes_Is_to_Know_How_Big_They_Really_Are)
