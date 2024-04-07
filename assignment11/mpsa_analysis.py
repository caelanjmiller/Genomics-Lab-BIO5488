# Complete the code as indicated in the comments.

# I strongly recommend working out individual sections of code
# in an interactive format like jupyter notebook. Then paste your
# code into the appropriate sections below and run this script to generate the necessary # output.

# needed imports
import pandas as pd

# read in the data
df = pd.read_csv("mpsa_data.csv", sep=",")

# visualize the dataframe structure
# Take note of the columns.
print("df.head():\n", df.head(), "\n")


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
# Use this list to track your decile counts
decile_counts = [0] * 10

# Extract the list of PSI values from the dataframe:

### YOUR CODE HERE

# Iterate through the list of PSI scores. For each
# score, determine its decile then increment the decile count.
# You can use get_decile_index to obtain the correct
# decile to increment in decile_counts.

### YOUR CODE HERE

# Print your output
print("Use the output below to answer question 1.")
print("Decile distribution:\n", decile_counts)

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


# Create list of GC deciles:
gc_decile_counts = [0] * 10

### YOUR CODE HERE

# Print your output:
print("Use the output below to answer question 2.")
print("Decile distribution:\n", gc_decile_counts)


#### Task 4: Print out GC sequences that are in the top decile ####
# Hint: Use print() to show either the relevant rows of the dataframe
# or create a list of sequences in the top decile and print the list.

### YOUR CODE HERE

#### Task 5: Check the PSI of pathogenic variants ####

# Read in list of pathogenic sequences
file = open("brca2_mutations.txt")
brca2_muts = []
for line in file:
    brca2_muts.append(line.strip("\n"))
file.close()

# Extract PSI scores for the pathogenic sequences
# Hint: the pandas function .isin() will be helpful here.
# Alternately, you can use a for loop to create a list of
# the correct psi scores.

### YOUR CODE HERE

# Calculate the distribution among deciles
# Similar to Tasks 2 and 3.
brca2_decile_counts = [0] * 10  # list to hold counts

# Again use get_decile_index for each score
# Use that index to increment the decile count

### YOUR CODE HERE


# Print out the distribution
print("Use the output below to answer question 4.")
print("Decile distribution:\n", brca2_decile_counts)

# Calculate the mean PSI for the pathogenic sequences and print
# the result:

### YOUR CODE HERE
