Please provide the exact command line arguments you used to generate your results.
{Your command with specific file names for finding_gaps.py}
# Find Sequence Length of each FASTA
python3 finding_gaps.py hg38-chr22.fa length
python3 finding_gaps.py t2t-chr22.fa length
# Nucleotide Frequencies
python3 finding_gaps.py hg38-chr22.fa frequency
python3 finding_gaps.py t2t-chr22.fa frequency 
# Dinucleotide Frequencies
python3 finding_gaps.py hg38-chr22.fa difrequency
python3 finding_gaps.py t2t-chr22.fa difrequency
-
{your command with specific file names for building index files and aligning with bowtie2}
bowtie2-build t2t-chr22.fa t2t-chr22-index
bowtie2 --threads 16 -x t2t-chr22-index -U test-500k.fq -S t2t_chr22_alignment.sam 2> t2t_chr22_alignment_report.txt
bowtie2-build hg38-chr22.fa hg38-chr22-index
bowtie2 --threads 16 -x chr22_index -U test-500k.fq -S hg38_chr22_alignment.sam 2> hg38_chr22_alignment_report.txt
-
{your commands with specific file names for identifying duplicate reads with samtools}

-
Part 1:
Question 1:
{size of chromosome 22 in hg38 and CHM13}
HG38: 50818468 bp
HCM13 (T2T): 51324926 bp
-
Question 2:
{nucleotide frequencies for hg38 and CHM13}
HG38
A:0.265
C:0.234
G:0.236
T:0.265

CHM13
A:0.272
C:0.229
G:0.230
T:0.270
-
{dinucleotide frequencies for hg38 and CHM13}
HG38
AA:0.079
AC:0.051
AG:0.075
AT:0.060
CA:0.076
CC:0.067
CG:0.016
CT:0.074
GA:0.062
GC:0.055
GG:0.068
GT:0.051
TA:0.047
TC:0.061
TG:0.077
TT:0.080

CHM13
AA:0.082
AC:0.051
AG:0.070
AT:0.068
CA:0.075
CC:0.064
CG:0.016
CT:0.073
GA:0.063
GC:0.050
GG:0.065
GT:0.051
TA:0.051
TC:0.064
TG:0.077
TT:0.077
-
{explanation}

-
Question 3:
{number of gaps in hg38 and CHM13}

-
{please list all gap lengths from hg38 and CHM13 that you found here}

-
Part 2:
Question 4:
{number of mappable reads for hg38 and CHM13}
HG38 - 
CHM13 - 
-
{number of non-mappable reads for hg38 and CHM13}
HG38 - 
CHM13 - 
-
{number of uniquely-mapped reads for hg38 and CHM13}

-
{exploration of similarities/differences}

-
Question 5:
{percent of uniquely-mapped reads for hg38 and CHM13}

-
{percent of multi-mapped reads for hg38 and CHM13}

-
{explanation}

-
Question 6:
{fraction of de-duplicated reads in total mappable reads for hg38 and CHM13}

-
{comparison}

-
Part 3:
{time for initial BLAST job}

-
Question 1:
{Explanation}

-
Question 2:
{Number of hits}

-
Question 3:
{Genus species, 
Score, % identity}

-
Question 4:
{Number of hits, explanation}

-
Question 5:
{Genus species}

-
Question 6:
{Score, % identity}

-
Question 7:
{Number of hits, explanation}

-
Question 8:
{Explanation}

-
Question 9:
{Answer, Explanation}

-
Question 10:
{}

-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}

-



