Usage:
{How to run neutral_rate.py}
neutral_rate.py PRE1.aln PRE1.fa
-
Question 1:
{Fraction of wobble positions conserved}
0.4623
-
Question 2
{Cutoff for 10bp sequence conserved}
7 bp
{Explanation}
I utilized the percent point function (ppf) for the binomial class, where I defined .4623 as the probability
of neutral conservation, 10 bp as the n, and .95 as 
num_bp = int(scipy.stats.binom.ppf(.95, 10, .4623))
-
Question 3:
{Number of regions}
{S_cer_consreved.txt output}
-
Question 4:
{Comparing conserved regions to full promoter region}
-
Question 5:
{Do binding sites point to proteasome assembly}
-