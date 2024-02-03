Question 1
{Your command with specific file names}
Bowtie Alignment:
bowtie2 --threads 16 -x chr22_index -U reads.fq -S chr22_alignment.SAM 2> alignment_report.txt

FASTQ Parser:
python3 fastq_parser.py <FASTQ> <PRINTOUT>
- Nucleotide Frequency 
python3 fastq_parser.py reads.fq frequency > reads_frequency.txt
- Dinucleotide Frequency
python3 fastq_parser.py reads.fq difrequency > reads_difrequency.txt

Nucleotide Count:
python3 nucleotide_count.py <FASTA> <CUSTOM FASTA HEADER> <PRINTOUT>
- chr22 Nucleotide Frequency
python3 nucleotide_count.py chr22.fa chr22 frequency > chr22_frequency.txt
- chr22 Dinucleotide Frequency
python3 nucleotide_count.py chr22.fa chr22 difrequency > chr22_difrequency.txt

Nucleotide Enrichment:
python3 nucleotide_enrichment.py <FILE 1> <FILE 2> 
- Nucleotide Enrichment
python3 nucleotide_enrichment.py reads_frequency.txt chr22_frequency.txt
- Dinucleotide Enrichment
python3 nucleotide_enrichment.py reads_difrequency.txt chr22_difrequency.txt

{Number of uniquely mapped reads}
7115 (26.69%)

{Number of multi mapped reads}
10126 (37.98%)

{Number of unmapped reads}
9421 (35.33%)

-
Question 2
{What is enriched in this dataset}
{Single nucleotide enrichment scores}
A:1.034
C:1.086
G:1.114
T:0.785

{Dinucleotide enrichment scores}
AA:1.165
AC:0.941
AG:0.733
AT:1.328
CA:0.895
CC:1.045
CG:3.938
CT:0.703
GA:1.302
GC:1.148
GG:1.191
GT:0.725
TA:0.702
TC:1.18
TG:0.831
TT:0.481

In this dataset, the nucleotides cytosine, guanine, & adenine are most enriched (fold greater than 1); 
this is also reflected in the enriched dinucleotides (e.g. GC, CG, CC, etc.). Of these dinucleotides,
CG (at ~4 fold greater in our sequencing reads) is the most enriched, followed by AA, CC, AT, GA, TC, GC, & GG (~1 fold change).

{Description of how enrichment scores were calculated}
These enrichment scores were produced by first calculating the individual nucleotide and dinucleotide
frequencies of both reads.fq & chr22.fa. I then took these scores and divided the sequencing reads frequencies by the frequencies 
for chr22. 


{Assay, Explanation}
Based upon the GC enrichment in our dataset, it is likely that this data is a result of a ChIP sequencing experiment. 
The areas of interest in a ChIP sequencing experiment (promoter regions, enhancers, TF binding sites, etc.) are quite GC-rich, which would explain the
CG nucleotide and dinucleotide enrichment seen in our dataset. Additionally, per Bowtie2's alignment report, many reads were mapped multiple times (multi mapped) and this 
could be due to the shared binding sites across enhancers, promoter regions, etc leading to sequencing reads being mapped across the reference multiple times.
-
Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
-
