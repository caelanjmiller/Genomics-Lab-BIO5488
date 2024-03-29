from sys import argv, exit
import numpy as np
from scipy import stats
import math

A: int = int(argv[1])
B: int = int(argv[2])
C: int = int(argv[3])
D: int = int(argv[4])


"""
Python script to calculate a 2-Tailed Fisher's exact test, returning:
1. Odds Ratio
2. P-Value
3. Upper Limit of 95% CI
4. Lower Limit of 95% CI

Usage: python3 fisher_test.py A B C D 
A: Number of 1 calls for alternate allele
B: Alternate allele count - gnomAD
C: Number of 0 calls for reference allele
D: Total allele number on gnomAD - alternate allele count (B)
"""

contingency_table = np.array([[A, B], [C, D]])
results = stats.fisher_exact(contingency_table)
upper_limit_CI = round(
    math.e ** (np.log(results[0]) + 1.96 * (math.sqrt(1 / A + 1 / B + 1 / C + 1 / D))),
    4,
)
lower_limit_CI = round(
    math.e ** (np.log(results[0]) - 1.96 * (math.sqrt(1 / A + 1 / B + 1 / C + 1 / D))),
    4,
)


print(f"Odds Ratio: {round(results[0], 4)}")
print(f"P-value: {round(results[1], 4)}")
print(f"Upper Limit of 95% CI: {upper_limit_CI}")
print(f"Lower Limit of 95% CI: {lower_limit_CI}")

if len(argv) != 5:
    print(__doc__)
    exit(1)
