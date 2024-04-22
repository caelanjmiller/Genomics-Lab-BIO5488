Part I


Question 1:
{Exact command to run count_gv.py}
python3 gv_count.py SNV_indel.biallelic.vcf sv.reclassed.filtered.vcf
{Output count table}
type    count
DEL     1664
DUP     530
BND     2690
MEI     1200
INV     95
INDEL   479363
SNV     4086192
-
Question 2:
{Proportion of genomic variants are SVs}
0.0014 - Calculated as sum(SVs) / total count of SVs & SNVs & indels
-
Question 3:
{Describe the spectrum of SVs for the individual NA12878}
The spectrum of large genomic structural variations in NA12878 is diverse, with the majority existing as 
breakend variations, followed by large deletions & mobile element insertions. 
-
Question 4:
{Place the three histograms in your submission directory}

{Describe the distributions observed in each of the histograms}
The distribution for INDELs tends to have a positive skew (skewing towards 0bp in length), with the majority
of INDELs being 50 bp or less in size. The distribution for MEIs also has a positive skew (skewing towards 0 bp), but has an overall
longer tail that include lengths up to 6000 bp. The distribution for larger deletions was quite varied as the spread between datapoints was 
quite large, some deletions could be as short as 70 bp while others could be up to 50,000 bp long. 
{Speculate how the length distribution might differ if we limit the data to exonic indels?}
I would imagine that if you limited the data to just exonic indels, the lengths would drastically shrink in overall magnitude as the evolutionary tolerance
for this potential lack of sequence would be quite low due to potential deleterious effects, which would in turn reduce fitness and
lead to a lower frequency within a population.
-
Part II

Question 5:
{Exact command to run quantify_genotype.py}
python3 quantify_genotype.py SNV_indel.biallelic.vcf
{Output count table}
type    count
homozygous_ref  3518113
heterozygous    2922423
homozygous_alt  1643132
missing 95900
-
Question 6:
{Does the difference in the number of homozygous alternate (or non-reference homozygous) and heterozygous SNVs and indels, make biological sense? Why, or why not}
Yes, they do make some biological sense. Homozygous reference alleles are present in the greatest frequency followed by heterozygous
alleles and this makes sense given that if we follow a simple dominance model (for simplicity), that genes that confer a fitness value would be 
passed along within a given population due to natural selection. Genes that would be otherwise be homozygous recessive (e.g homozygous alternate) and confer
a disadvantage would be present at a lower frequency compared to their wild type homozygous and heterozygous counterparts.
-
Question 7:
{How many variants clearly violate the rules of Mendelian segregation? (Autosome only)}

-
Question 8:
{Describe four potential reasons that could explain the Mendelian violations.}

-
Question 9:
{How many variants now violate mendelian segregation after filtering? (Autosome only)}
-
Comments:
{Things that went wrong or you can not figure out}
- Parsing the VCF took sometime to figure out due to inconsistency in fields
