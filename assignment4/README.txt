Part 1.0
{Command for running analyze_WGBS_methylation.py}
python3 analyze_WGBS_methylation.py BGM_WGBS.bed 
{Copy your output files BGM_WGBS_CpG_methylation.bed, BGM_WGBS_methylation_distribution.png, and BGM_WGBS_CpG_coverage_distribution.png to your submissions directory}
-
Question 1:
{What does DNA methylation look like across chromosome 21?}
DNA methylation in chromosome 21 tends to skew right (towards average methylation of 0), but there is a significant peak towards 
an average methylation level of 1.0
Question 2:
{What does the CpG coverage look like across chromosome 21?}
The data tends to skew right (towards a read coverage of 0), indicating a lot of the chromosome does not have a lot of coverage in 
terms of sequence space, this is further complicated by the fact that we also likely did not have great depth in those reads so we may not
be able to be fully confident on the identity of the called base (in our case since we're focused on bisulfite sequencing, C & T)
-
Question 2.1:
{What fraction of the CpGs have 0X coverage?}
0.084
-
Part 1.1
{Command for creating a bed file with the average CpG methylation level in each CGI.}
# Create bed file for methylated CpGs that overlap with CGI intervals 
bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b CGI.bed -wb | sort -k 1,1 -k2,2n > CpG_CGI_overlap.bed
# Create bed file with average methylation in each CGI
bedtools groupby -i CpG_CGI_overlap.bed -g 5-8 -c 4 -o mean > WGBS_CGI_methylation.bed
-
Part 1.2
{Command for plotting the distribution of average CGI methylation levels}
{Copy analyze_CGI_methylation.py and WGBS_CGI_methylation_distribution.png to your submissions directory}
python3 analyze_CGI_methylation.py WGBS_CGI_methylation.bed
-
Question 3:
{What does DNA methylation look like for CpGs in CGIs? How does it compare to all the CpGs on chromosome 21?}
DNA methylation in CGIs tends to follow a similar trend as seen with CpGs in all of chromosome 21
-
Part 1.3.0
Gene promoters
{Command for generating the promoter bed file}
python3 generate_promoters.py refGene.bed
{Justification for promoter definition}
So many sources [https://www.addgene.org/mol-bio-reference/promoters/, https://www.frontiersin.org/articles/10.3389/fbioe.2019.00305/full, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6848157/, amongst others] 
claim that promoter regions can range in length from 100 to 1000 base pairs. I chose 1000 bp because it is a long enough stretch to potentially encapsulate promoter elements
such as the actual core promoter along with proximal and distal elements. 
-
Promoter-CGI and non-promoter-CGI
{Commands for generating promoter-CGI and non-promoter-CGI bed files}
# Generating promoter-CGI
bedtools intersect -a WGBS_CGI_methylation.bed -b refGene_promoters.bed -wb | sort -k 1,1 -k2,2n > promoter_CGI.bed
# Generating non-promoter-CGI bed files
bedtools intersect -a WGBS_CGI_methylation.bed -b promoter_CGI.bed -v | sort -k 1,1 -k2,2n > non_promoter_CGI.bed
{Justification for overlapping criteria}
So for creating the files, I wanted to parse out gene coordinates that simultaneously existed in my CGI methylation data & regions defined
as gene promoters (for promoter_CGI.bed). For the non_promoter regions, I wanted regions that DID NOT overlap (hence the -v flag) between
my defined promoter regions and the CGI methylation data.  
{Commands for calculating the average CpG methylation for each promoter-CGI and non-promoter-CGI}
# Generating average promoter CGI methylation bed file
bedtools groupby -i promoter_CGI.bed -g 1-4 -c 5 -o mean | sort -k 1,1 -k2,2n > average_promoter_CGI_methylation.bed 
# Generating average non-promoter CGI methylation bed file
bedtools groupby -i non_promoter_CGI.bed -g 1-4 -c 5 -o mean | sort -k 1,1 -k2,2n > average_non_promoter_CGI_methylation.bed
{Commands for running analyze_CGI_methylation.py on average_promoter_CGI_methylation.bed and average_non_promoter_CGI_methylation.bed}
python3 analyze_CGI_methylation.py average_promoter_CGI_methylation.bed
python3 analyze_CGI_methylation.py average_non_promoter_CGI_methylation.bed
-
Question 4:
{How do the DNA methylation profiles of promoter-CGIs and non-promoter-CGIs differ?}
There appears to be drastic reduction in the methylation levels surrounding promoter regions (large number of average methylation levels equaling & skewing towards 0)
versus non-promoter regions (skewed towards the opposite direction; towards average methylation level of 1)
-
Part 1.3.1
{Commands for calculating CpG frequency for each promoter type}
python3 nucleotide_count.py promoter_CGI.fasta
python3 nucleotide_count.py non_promoter_CGI.fasta
{CpG frequencies for each promoter type}
Including CpG frequencies as separate files - lots of data:
Promoter CpG Frequency: 0.111
Non Promoter CpG Frequency: 0.089
-
Question 5:
{What is a possible biological explanation for the difference in CpG frequencies?  Interpret your results from parts 1.3.0 and 1.3.1: what are the “simple rules” for describing regulation by DNA methylation in promoters?}
GC rich regions are more thermostable and less prone to interruption in processes such as DNA transcription via RNA polymerase. Extrapolating
out to these promoters that are GC rich, this makes sense biologically as important areas where transcription needs to occur (aka your genes), one would
want these processes to occur with great efficiency. It has been demonstrated that these CpG islands surrounding promoter regions are 
often unmethylated and enriched for chromatin modifications leading to active transcription [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1345710/]. 
With these promoter regions being less prone to methylation (methylation can lead to eventual deamination of cytosine to thymine), despite
often being GC rich, this could lead one to infer that methylation may serve as a form of repression if a GC rich promoter is methylated. 
Methylation could impede DNA transcription by steric hindrance, with the bulky methyl groups blocking transcriptional machinery 
or these methyl groups could recruit other protein complexes and form the transcriptionally inactive, heterochromatin. 
-