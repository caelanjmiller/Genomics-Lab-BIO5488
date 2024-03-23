Due date: 03/29/2024 12 am (by end of Friday)


Part 1
A. Which GATK version is used in this dataset?
4.2.2.0

B1. How many samples are there in this VCF file?
500

B2. Please also write down the bcftools command you used to count the number of samples in this VCF file.
bcftools query -l example.vcf | wc -l 

C. Which reference genome is used in this VCF file?
GRCh38

Part 2
A1. What is the chromosome number and position of this variant?
chr17  28381741

A2.  Please write down the bcftools command you used to find the chromosome number and position in the README.txt file
bcftools query -f '%CHROM  %POS' example.vcf 

B. What’s the alternative allele of this variant?
G
bcftools query -f '%ALT' example.vcf

C. What’s the reference allele of this variant?
A
bcftools query -f '%REF' example.vcf

D1. What is the count of different genotypes present for this variant? Please print the genotypes out and calculate the number for each genotype. 
341  ./.
157  0/0
2  0/1

D2. Please also write down the bcftools command you used to count different genotypes present for this variant.
bcftools query -f '[ %GT\n]' example.vcf | sort | uniq -c

Part 3
A. In which gene does this genetic variant lie?
sterile alpha and TIR motif containing 1 (SARM 1)

B. How many alternative alleles are present in the gnomAD v3.1.2 database?
28

C. What’s the overall minor allele frequency in the gnomAD v3.1.2 database?
1.84e-4

D. What’s the variant type of this variant? (pLoF/missense/synonymous)
missense

 
E. Have any human diseases been reported to be associated with this gene's mutations? Do you believe that this variant is pathogenic to the disease? Please utilize information from ClinVar to make your judgment. Additionally, what is the known functional role of this gene according to the literature?
Congenital defect of folate absorption, Amyotrophic lateral sclerosis & general inborn genetic diseases have been identified as human diseases
associated with this gene's mutations. Based upon the nature of the mutations (missense), 

Part 4
A. Please write your Python script for the Two-tailed Fisher's exact test, and name it ‘fisher_test.py’.


B. Your Python script needs to be able to print out the following information on the screen: 
1) odds ratio, 
2) p-value,
3) upper limit of 95% confidence interval, 
4) lower limit of 95% confidence interval.




C. Copy and paste the calculated odds ratio, p-value, and upper and lower limits of 95% confidence interval here.




[END]