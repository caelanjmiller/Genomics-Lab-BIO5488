Bio5488 2024 Spring Assignment 11

Due: Friday April 12th 2024, by 11:59 pm

What to turn in:
1. mpsa_analysis.py 
2. README.txt
3. output.txt 

Please provide the exact command line arguments you used to generate your results.
{How to run mpsa_analysis.py}
python3 mpsa_analysis.py mpsa_data.csv brca2_mutations.txt
-
Question 1: Looking at the distribution of scores by PSI decile, into which two deciles to most PSI sequences fall? What does this result suggest about the efficiency 5’ splice sites?

Most PSI sequences fall into the 0 & 9th deciles; this would suggest that 5' splice sites are quite efficient, any mutation in the splice site obliterates any potential reporter read out. 
-
Question 2: Consider the number of ‘GC’ splice site sequences that fall into the top decile. Are ‘GC’ splice sites more or less likely than ‘GU’ splice sites to be spliced with 90-100% efficiency? What fraction of the splice sites in the top PSI decile are ‘GC’?

I'd say that 'GC' splice sites are more likely than 'GU' splice sites to be spliced with that high degree of efficiency (90-100%) as indicated in the data, due to being such a rare event. The fraction of GC non-canonical sites as a whole is ~50% (15414 rows of data), and of this, only 20 were in the top decile (9). So 20/30483 (the entire MPSA dataset), equals ≈.0006561.
-
Question 3: What fraction of the GC sequences in the top decile have G’s at both -1 and +5 positions? Sequences with G’s at -1 and +5 make up about 6% of all possible sequences. Are these sequences enriched in among sequences in top PSI decile?

Fraction of top decile GC sequences that have G's at both -1 & +5 positions: 0.2
These sequences are enriched in the GC sequences in the top PSI decile (0.2 >> 0.06)

-
Question 4: (A) Of 41 pathogenic mutations, what fraction have PSI scores in the top decile? In the entire MPSA dataset, 1.7% of sequences have scores in the top decile. Why would the mutation set differ in this fraction? (B) Among the wild-type versions of the pathogenic variants, the mean PSI score is > 99%. What is the mean PSI among the pathogenic variants?

A) Fraction of pathogenic mutations that have a PSI score in the top decile: 0.1591
This number would differ in the mutation set due to the fact that mutations in these splice sites (any intron) often lead to pathogenesis, so having an enrichment in this set of known pathogenic variants versus the ENTIRE dataset makes sense, hence the determination of pathogenicity for a given splice variant within the smaller dataset.

B) Mean PSI value of BRCA2 pathogenic variants: 100.0055

-