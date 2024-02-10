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
# Generation of Gap Length Histograms
python3 finding_gaps.py hg38_chr22.fa graph
python3 finding_gaps.py t2t_chr22.fa graph
-
{your command with specific file names for building index files and aligning with bowtie2}
bowtie2-build t2t-chr22.fa t2t-chr22-index
bowtie2 --threads 16 -x t2t-chr22-index -U test-500k.fq -S t2t_chr22_alignment.sam 2> t2t_chr22_alignment_report.txt
bowtie2-build hg38-chr22.fa hg38-chr22-index
bowtie2 --threads 16 -x chr22_index -U test-500k.fq -S hg38_chr22_alignment.sam 2> hg38_chr22_alignment_report.txt
-
{your commands with specific file names for identifying duplicate reads with samtools}
# Collates read information, fixes/repairs read information & sorts reads 
samtools collate -@ 4 -O -u hg38_chr22.sam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u -o hg38_chr22_positionsort.sam
samtools collate -@ 4 -O -u t2t_chr22.sam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u -o t2t_chr22_positionsort.sam

# Marks duplicates in sorted bam and outputs similar stats to MarkDuplicates (Picard) [https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard]
samtools markdup -@ 4 -f hg38_chr22_stats_sam.txt -S -d 2500 --mode s --include-fails hg38_chr22_positionsort.sam hg38_chr22_markdup.sam
samtools markdup -@ 4 -f t2t_chr22_stats_sam.txt -S -d 2500 --mode s --include-fails t2t_chr22_positionsort.sam t2t_chr22_markdup.sam
-
Part 1:
Question 1:
HG38: 50818468 bp
CHM13 (T2T): 51324926 bp
-
Question 2:
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
The nucleotide and dinucleotide frequencies of the two assemblies are relatively similar (this makes sense given they are both assemblies of the human genome!). However, 
where they mainly differ is the frequency of thymine (T), adenine (A) & guanine (G) and the respective dinucleotides that contain these nucleotides.
The CHM13 assembly (aka T2T [https://www.science.org/doi/10.1126/science.abj6987]) aimed to fill in the gaps seen in HG38 (represented by N), ultimately correcting many of the sequencing errors in the prior assembly and 
introducing resolved sequence space. The differences between the two assemblies as aforementioned with HG38 containing a lot of gaps (N), and since N is just
any nucleotide identity, any resolution on the actual identity of those gapped sequences will influence the frequencies. On an additional note, with the new sequence resolution near
telomeric regions (which are known by a consensus sequence of (TTAGGG)n [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC120798/]), this may be the reason
why one would see an increase/difference in the frequencies of three aforementioned nucleotides (A, G, T)
-
Question 3:
{number of gaps in hg38 and CHM13}
HG38: 49 gaps
CHM13: 0 gaps
-
{please list all gap lengths from hg38 and CHM13 that you found here}
CHM13 aka T2T: 0 gaps
HG38:
Length: Count
10510000:1
50000:18
100:15
7856:1
2477:1
100000:2
23171:1
1131:1
557:1
1:5
494:1
1500:1
10000:1
-
Part 2:
Question 4:
{number of mappable reads for hg38 and CHM13}
HG38:
500000 reads; of these:
  500000 (100.00%) were unpaired; of these:
    133832 (26.77%) aligned 0 times
    309685 (61.94%) aligned exactly 1 time
    56483 (11.30%) aligned >1 times
73.23% overall alignment rate

CHM13:
500000 reads; of these:
  500000 (100.00%) were unpaired; of these:
    101096 (20.22%) aligned 0 times
    316543 (63.31%) aligned exactly 1 time
    82361 (16.47%) aligned >1 times
79.78% overall alignment rate
-
{exploration of similarities/differences}
For CHM13, I had a better overall alignment of my sequencing reads to the assembly than HG38. This is likely due to the 
aforementioned sequence gap filling that occurred with this assembly so the actual sequence space that could be mapped back with 
reads is significantly larger (this is most apparent in the nearly 7% reduction in nonmapped reads (aligned 0 times)). 
-
Question 5:
{percent of uniquely-mapped reads for hg38 and CHM13}
[T2T] 316543 / 500000 = .6331
[HG38] 309685 / 500000 = .6194
-
{percent of multi-mapped reads for hg38 and CHM13}
[T2T] 82361 / 500000 = .1647
[HG38] 56483 / 500000 = .1130
-
{explanation}
I think the explanation for the increased percentage of both uniquely and multimapped reads are due to the gapless nature of 
the T2T assembly. Having defined sequences allows bowtie2 to map the reads back onto the reference (the respective assembly) with
greater confidence and accuracy and from the publication about T2T (linked above), quite a few coding regions were potentially elucidated
with the new assembly and it's possible they share some similarity in already established (sequenced) regions, leading to more multi mapped reads.
The addition of new coding regions as well elucidates new genomic regions, allowing for a greater increase of uniquely mapped reads as well
-
Question 6:
{fraction of de-duplicated reads in total mappable reads for hg38 and CHM13}
-T2T-
COMMAND: samtools markdup -@ 4 -f t2t_chr22_stats_sam.txt -S -d 2500 --mode s --include-fails t2t_chr22_positionsort.sam t2t_chr22_markdup.sam
READ: 500000
WRITTEN: 500000
EXCLUDED: 101096
EXAMINED: 398904
PAIRED: 4
SINGLE: 398900
DUPLICATE PAIR: 0
DUPLICATE SINGLE: 96059
DUPLICATE PAIR OPTICAL: 0
DUPLICATE SINGLE OPTICAL: 67
DUPLICATE NON PRIMARY: 2
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 96059
DUPLICATE TOTAL: 96061
ESTIMATED_LIBRARY_SIZE: 0

-HG38-
COMMAND: samtools markdup -@ 4 -f hg38_chr22_stats_sam.txt -S -d 2500 --mode s --include-fails hg38_chr22_positionsort.sam hg38_chr22_markdup.sam
READ: 500000
WRITTEN: 500000
EXCLUDED: 133832
EXAMINED: 366168
PAIRED: 8
SINGLE: 366160
DUPLICATE PAIR: 2
DUPLICATE SINGLE: 90980
DUPLICATE PAIR OPTICAL: 2
DUPLICATE SINGLE OPTICAL: 65
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 90982
DUPLICATE TOTAL: 90982
ESTIMATED_LIBRARY_SIZE: 0

Looking at the stat read-outs for samtools markdup, I would calculate the % of deduplicated reads out of total mappable reads as:
[T2T]  398900 (total number of deduplicated single reads) / 500000 (total amount of reads) = .7978
[HG38] 366160 (total number of deduplicated single reads) / 500000 (total amount of reads )= .7323
-
{comparison}
Again, as I have noticed with the bowtie2 alignment even prior to the deduplication processing via samtools, I 
can tell that I obtained more overall mapped reads for the T2T assembly. On the sequencing preparation side, it could be likely as well
that there is variation with library generation and this could lead to having overall bad (duplicated reads, sequencing errors) read alignments
back onto the reference genome assembly
-
Part 3:
~ 2 minutes

-
Question 7:
Utilizing the non-redundant (nr) database ensures that you are not getting multiple hits for the same protein many different times 
in the same submitted organisms. The database has been set up in such a manner that each submitted reference 
genome/proteome has been submitted once and is of the lowest taxonomic classification possible. This is likely to prevent database biasing towards
the number of submitted reference genomes/isolates (e.g. E.coli K12 has a very large number of submitted assemblies but only 1 is 
a reference utilized in BLAST results, unless it is a different strain/isolate).

-
Question 8:
I got ~412 putative hits on 500 subject sequences with the BLOSUM62 Scoring Matrix

-
Question 9:
{Genus species, 
Score, % identity}
Torulaspora globosa
871 - Max & Total Score
60.44% Percent Identity

-
Question 10:
I set my Max target sequences to 500 and received ~417 putative hits with an e-value < 1 with the BLOSUM80 Scoring Matrix; I got more hits than with the BLOSUM62
In theory, I should have gotten less hits, given that BLOSUM80 assumes more closely related alignments and thus would penalize
for sequences that do not align (use of this matrix would assume prior knowledge/inference of evolutionary results [https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Use_scoring_matrices.html])

-
Question 11:
Torulaspora globosa
902 - Max & Total Score
60.30% Percent Identity

-
Question 12:
698 & 806 - Max & Total Score
74.10% Percent Identity
-
Question 13:
Setting my Max target sequence to 500 and lowering the gap costs to 7, I received ~396 hits on 321 subject sequences with an e-value < 1 with the BLOSUM62
Scoring Matrix (compared to ~412 hits on 346 subject sequences). In theory, I believed I should have gotten more hits as I am reducing the overall
gap penalty score, which coupled with a scoring matrix that assumes distantly related proteins (BLOSUM62), would return more hits as 
the criteria for an optimal alignment (according to the algorithm parameters rather than being biologically relevant) has been reduced
-
Question 14:
Changing the gap penalty parameters allows for orthologs (especially when using the BLOSUM62 scoring matrix - used for distant alignments), that have higher gaps to
potentially obtain higher overall scores (because of lower penalty against their overall alignment score) 
and thus pop up higher on the returned subject list

-
Question 15:
I would expect the search time to be longer as the search becomes more sensitive (shorter sequence substring to be searched) but the 
actual search results that are returned are less specific (and potentially less biologically relevant)

-
Question 16:
Yes, but this depends on query size and various scoring parameters but it is very slow compared to potentially running 
the BLAST+ binaries (with the various databases downloaded) on a dedicated server

-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}

-



