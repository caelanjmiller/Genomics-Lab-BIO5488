#!/usr/bin/env python3


"""
Python script to analyze MPSA data for BRCA2
Usage: python3 mpsa_analysis.py <MPSA DATA CSV> <BRCA VARIANTS TXT>
"""

import pandas as pd
import numpy as np
from sys import argv
from pathlib import Path
from collections import Counter

MPSA_CSV: Path = Path(argv[1])
BRCA_TXT: Path = Path(argv[2])


if len(argv) != 3:
    print(__doc__)
    exit(1)

# read in the data
MPSA_DF = pd.read_csv(MPSA_CSV, sep=",")


# A function to determine PSI decile:
def get_decile_index(psi, minumum=0, bin_width=10):
    """Returns the decile  for measurements scaled 0 to 100. PSI of
    >= 100 is assigned the last decile."""
    index = int((psi - minumum) // bin_width)
    if index > 9:
        index = 9
    return index


#### TASK 1: Write a helper function ####
# Complete the function below. Recall that the
# MPSA sequences are 9 nt long that begin at the -3 position.
# They are represented in python as strings that begin
# with string position 0. See fig. 2 in the instructions.

# The function should take a single splice site position
# as input (e.g., -2) and return the python string coordinate.


def splice_site_coordinate_conversion(splice_site_index: int) -> int:
    """Convert splice site coordinates into 0-based indexing"""
    if splice_site_index < -3 or splice_site_index > 6:
        return None
    elif splice_site_index == -3:
        return 0
    else:
        return splice_site_index + 2


#### TASK 2: Examining the distribution of PSI values ####

# How are PSI values distributed? Complete the code to calculate the number of
# 5' ss sequences that fall into each PSI decile, e.g., PSI values of 0-10%
# fall in the first decile, etc.

# You will need to (1) Check each sequence and determine which decile
# its PSI falls in, then (2) increment the sequence count for that decile.

# To help you, you may use the get_decile_index() function to covert a PSI score
# to a decile index. The function takes a PSI as input and returns a value from 0 to 9.
# 0 indicates the 1st decile, 1 is the 2nd decile, etc.
# We start at 0 so that the value can be used as an index for the
# list of decile counts.

# Use the code and hints below as a starting point or choose your own approach.

# Initialize a list of counts for 10 deciles.
# Values in the list represent number of sequences
# with a PSI within [0:10), [10:20), up through [90:100]
# Due to experimental noise some sequences have a PSI > 100.
# These are assumed to be equivalent to PSI 100.

# Extract the list of PSI values from the dataframe:

psi_values: list = MPSA_DF["psi"].to_list()
decile_indices: list = [get_decile_index(psi) for psi in psi_values]
decile_indices_counter: dict = dict(Counter(decile_indices))

# Print your output
print("Use the output below to answer Question 1")
print(f"Decile indices: {dict(sorted(decile_indices_counter.items()))}")

#### TASK3: Examine the distribution of PSI values of GC splice sites.

# Extract splice sites with GC in +1,+2 position, get their PSI scores,
# then determine their distribution among PSI deciles.

# Hint 1: To identify GC sequences, use ss_to_string() to convert
# the +1 and +2 positions to python string coordinates. Then look for matches
# to a substring like this: my_string[pos1:pos2+1].

# Hint 2: Collect GC sequences and their PSIs either by
# creating a list of sequences and a list of PSIs from the
# dataframe df. Or you can create a new GC-specific dataframe
# by selecting only those rows of df with GC splice sites. This can
# done with a function like df.apply().

# Hint 3: To calculate the distribution among
# deciles, you can recyle your code from task 2 with a few small tweaks.

mpsa_data_dict: dict = dict(zip(MPSA_DF["splice_site"], MPSA_DF["psi"]))
gc_splice_site_psi: dict = {}
for splice_site, psi in mpsa_data_dict.items():
    position_one: int = splice_site_coordinate_conversion(1)
    position_two: int = splice_site_coordinate_conversion(2)
    if splice_site[position_one : position_two + 1] == "GC":
        gc_splice_site_psi[splice_site] = psi
    else:
        continue
gc_splice_site_psi_values: list = list(gc_splice_site_psi.values())
# Create list of GC deciles:
gc_splice_site_decile_indices: list = [
    get_decile_index(psi) for psi in gc_splice_site_psi_values
]
gc_splice_decile_indices_counter: dict = dict(Counter(gc_splice_site_decile_indices))

# Print your output:
print("Use the output below to answer question 2.")
print(f"Decile distribution: {dict(sorted(gc_splice_decile_indices_counter.items()))}")


#### Task 4: Print out GC sequences that are in the top decile ####
# Hint: Use print() to show either the relevant rows of the dataframe
# or create a list of sequences in the top decile and print the list.

top_decile_gc_splice_site_sequences: list = []
for splice_site, psi in gc_splice_site_psi.items():
    decile_index: int = get_decile_index(psi)
    if decile_index == 9:
        top_decile_gc_splice_site_sequences.append(splice_site)
    else:
        continue
print("GC Sequences in Top Decile:")
print(top_decile_gc_splice_site_sequences)

# Determine fraction of GC sequences in top decile that have G's at both -1 & +5 positions
gc_have_gg_counter: int = 0
for sequence in top_decile_gc_splice_site_sequences:
    position_one: int = splice_site_coordinate_conversion(-1)
    position_two: int = splice_site_coordinate_conversion(+5)
    if sequence[position_one] == "G" and sequence[position_two] == "G":
        gc_have_gg_counter += 1
    else:
        continue
print(
    f"Fraction of top decile GC sequences that have G's at both -1 & +5 positions: {round(float(gc_have_gg_counter / len(top_decile_gc_splice_site_sequences)), 4)}"
)

#### Task 5: Check the PSI of pathogenic variants ####

# Read in list of pathogenic sequences
brca2_mutations: list = [
    mutation.strip() for mutation in open(Path("brca2_mutations.txt")).readlines()
]

# Extract PSI scores for the pathogenic sequences
# Hint: the pandas function .isin() will be helpful here.
# Alternately, you can use a for loop to create a list of
# the correct psi scores.

# Calculate the distribution among deciles
# Again use get_decile_index for each score
# Use that index to increment the decile count


brca2_filtered_mutations_df = MPSA_DF[MPSA_DF["splice_site"].isin(brca2_mutations)]
brca2_psi_values = brca2_filtered_mutations_df["psi"].to_list()
brca2_decile_indices: list = [get_decile_index(psi) for psi in brca2_psi_values]
brca2_decile_indices_counter: dict = dict(Counter(brca2_decile_indices))

# Print out the distribution
print("Use the output below to answer question 4.")
print(
    f"BRCA2 Decile distribution: {dict(sorted(brca2_decile_indices_counter.items()))}"
)

# Calculate the mean PSI for the pathogenic sequences and print
# the result:

brca2_top_decile_psi_values: list = []
for psi in brca2_psi_values:
    decile_index: int = get_decile_index(psi)
    if decile_index == 9:
        brca2_top_decile_psi_values.append(psi)
    else:
        continue
brca2_mean_psi: float = round(np.mean(brca2_top_decile_psi_values), 4)
print(f"Mean PSI value of BRCA2 pathogenic variants: {brca2_mean_psi}")

print(
    f"Fraction of pathogenic mutations that have a PSI score in the top decile: {round(float(len(brca2_top_decile_psi_values) / len(brca2_mutations)), 4)}"
)
