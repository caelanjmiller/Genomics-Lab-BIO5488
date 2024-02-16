Part 1.0
{Command for running analyze_WGBS_methylation.py}
python3 analyze_WGBS_methylation.py BGM_WGBS.bed 
{Copy your output files BGM_WGBS_CpG_methylation.bed, BGM_WGBS_methylation_distribution.png, and BGM_WGBS_CpG_coverage_distribution.png to your submissions directory}
-
Question 1:
{What does DNA methylation look like across chromosome 21?}
-
Question 2:
{What does the CpG coverage look like across chromosome 21?}
-
Question 2.1:
{What fraction of the CpGs have 0X coverage?}
0.084
-
Part 1.1
{Command for creating a bed file with the average CpG methylation level in each CGI.}
# Create bed file for methylated CpGs that overlap with CGI intervals 
bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b CGI.bed -wb -sorted > CpG_CGI_overlap.bed
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
-
Part 1.3.0
Gene promoters
{Command for generating the promoter bed file}
python3 generate_promoters.py refGene.bed
{Justification for promoter definition}
So many sources [https://www.addgene.org/mol-bio-reference/promoters/, https://www.frontiersin.org/articles/10.3389/fbioe.2019.00305/full, amongst others] claim that promoter regions
can range in length from 100 to 1000 base pairs. I chose 500 bp because it is the middle ground between these two lengths; 100 bp may be too short for some long genes (multiple exons spanning kbp of sequence space) 
but 1000 bp may be too long for promoters in shorter genes. This approximation that I have taken is not comprehensive but aims to minimize losing out on either end of the spectrum in terms of capturing promoter
regions.
{Copy generate_promoters.py and refGene_promoters.bed to your submissions directory}
-
Promoter-CGI and non-promoter-CGI
{Commands for generating promoter-CGI and non-promoter-CGI bed files}
# Generating promoter-CGI
bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b refGene_promoters.bed -wb -sorted > promoter_CGI.bed
# Generating non-promoter-CGI bed files
bedtools intersect -a refGene.bed -b promoter_CGI.bed -v -sorted > non_promoter_CGI.bed
{Justification for overlapping criteria}
{Commands for calculating the average CpG methylation for each promoter-CGI and non-promoter-CGI}
{Commands for running analyze_CGI_methylation.py on average_promoter_CGI_methylation.bed and average_non_promoter_CGI_methylation.bed}
{Copy refGene_promoters.bed, promoter_CGI.bed, non_promoter_CGI.bed average_promoter_CGI_methylation.bed, average_non_promoter_CGI_methylation.bed, average_promoter_CGI_methylation.png and average_non_promoter_CGI_methylation.png to your submissions directory}
-
Question 4:
{How do the DNA methylation profiles of promoter-CGIs and non-promoter-CGIs differ?}
-
Part 1.3.1
{Commands for calculating CpG frequency for each promoter type}
{CpG frequencies for each promoter type}
-
Question 5:
{What is a possible biological explanation for the difference in CpG frequencies?  Interpret your results from parts 1.3.0 and 1.3.1: what are the “simple rules” for describing regulation by DNA methylation in promoters?}
-