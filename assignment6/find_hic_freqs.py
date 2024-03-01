#!/usr/bin/env python3
"""
TODO: write a doc string
<Your description of this script>
Usage: python3 find_hic_freq.py GENOMIC_BIN_MATRIX GM_MATRIX K562_MATRIX TAD_CSV
"""

# Loading recommended packages
from sys import argv, exit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

GENOMIC_BIN_MATRIX = Path(argv[1])
GM_MATRIX = Path(argv[2])
K562_MATRIX = Path(argv[3])
TAD_CSV = Path(argv[4])

# TODO: check the correct number of arguments were provide and exit the script if they were not
if len(argv) != 5:
    print(__doc__)
    exit(1)


# TODO: Read in the normalized matrix for the GM12878 and K562 cells
def parse_matrices(TSV: Path) -> np.ndarray:
    """Parse in TSV files as Numpy matrices"""
    numpy_matrix = np.loadtxt(TSV, delimiter="\t", dtype=float)
    return numpy_matrix


def parse_CSV(CSV: Path) -> pd.DataFrame:
    """Parse in CSV as pandas Dataframe"""
    dataframe = pd.read_csv(CSV, header=0)
    return dataframe


genomic_bin_matrix = parse_matrices(GENOMIC_BIN_MATRIX)
gm_matrix = parse_matrices(GM_MATRIX)
k562_matrix = parse_matrices(K562_MATRIX)
tad_dataframe = parse_CSV(TAD_CSV)


# TODO: Read in the list of genomic loci (genomic_bins.tsv) and the positional information of the TAD boundaries.
def calculate_genomic_bins(
    genomic_bin_matrix: np.ndarray, tad_dataframe: pd.DataFrame
) -> tuple:
    """Assign indices of genomic bins according to the range of TAD regions"""
    # Convert Dataframe to dict of lists
    tad_list: list = tad_dataframe.to_dict(orient="records")
    # Create dict that has tuple of ranges
    tad_dict: dict = {}
    for row in tad_list:
        row_values = list(row.values())
        tad_dict[row_values[0]] = tuple(row_values[1:])
    tad1_bins: list = []
    tad2_bins: list = []
    # Check to see if each genomic bin is in range and if so assign to appropriate TAD
    for index, genomic_bin in enumerate(genomic_bin_matrix):
        tad1_intervals: tuple = tad_dict[list(tad_dict.keys())[0]]
        tad2_intervals: tuple = tad_dict[list(tad_dict.keys())[1]]
        if genomic_bin >= tad1_intervals[0] and genomic_bin <= tad1_intervals[1]:
            tad1_bins.append(index)
        elif genomic_bin >= tad2_intervals[0] and genomic_bin <= tad2_intervals[1]:
            tad2_bins.append(index)
        else:
            continue
    # Create appropriate dictionary for data indices
    tad1_genomic_bin_data_indices: dict = {}
    tad2_genomic_bin_data_indices: dict = {}

    tad1_genomic_bin_data_indices[list(tad_dict.keys())[0]] = tuple(
        (min(tad1_bins), max(tad1_bins))
    )
    tad2_genomic_bin_data_indices[list(tad_dict.keys())[1]] = tuple(
        (min(tad2_bins), max(tad2_bins))
    )
    return tuple((tad1_genomic_bin_data_indices, tad2_genomic_bin_data_indices))


tad_genomic_bin_data_indices: tuple = calculate_genomic_bins(
    genomic_bin_matrix, tad_dataframe
)

# TODO: subset the genomic location of each TAD

# TODO: find the indices of the start and end locations of the TAD boundaries within the list of genomic locations.


## Part 1

# Complete the function plot_hic_map(). This function should take as input the matrix to plot, the TAD boundary indices to plot, the title, the file name, and the percentile cutoff for the maximum intensity (default to 95th percentile)

# TODO: write a function plot_hic_map() that will plot a heatmap of the contact matrices.
# As part of this function you should define the maximum intensity of the heatmap based on the n'th percentile of the data (do not harcode this).
# Identify the value from the input matrix that is the n'th percentile (n here is one of your inputs to the function)
# Plot a heatmap of the matrix (use a colormap in this plot. one way to fine control the colormap is the use the LinearSegmentedColormap.from_list() method to generate the colormap used in the plot).
# The minimum value should be 0 and the maximum should be set to the value of the nth percentile you generated.
# add rectangles to the heatmap where the TAD boundaries are if the TADs are present in the heatmap. If the length of the boundaries is larger than given bins, skip this process.
# add a title
# show the labels of x and y axis(i.e., the bin region)


def plot_hic_map(
    genomic_bin_data_indices: tuple,
    matrix_data: np.ndarray,
    title: str,
    filename: str,
    genomic_bin_matrix: np.ndarray,
    percentile: int = 95,
) -> None:
    """Create a heatmap of provided contact matrices"""
    # Unpack genomic bin locations for a given TAD
    tad1_genomic_data_index, tad2_genomic_data_index = genomic_bin_data_indices

    tad1_genomic_data_start, tad1_genomic_data_end = tad1_genomic_data_index["TAD1"]
    tad2_genomic_data_start, tad2_genomic_data_end = tad2_genomic_data_index["TAD2"]

    tad1_length = tad1_genomic_data_end - tad1_genomic_data_start
    tad2_length = tad2_genomic_data_end - tad2_genomic_data_start

    # Define genomic start and stop locations, and break up into 10 contiguous segements for x axis labeling
    labels_genomic_bins: list = (
        np.round(
            (np.linspace(min(genomic_bin_matrix), max(genomic_bin_matrix), 10) / 10**6),
            1,
        )
        .astype(str)
        .tolist()
    )
    final_labels: list = [f"{marker}MB" for marker in labels_genomic_bins]

    current_directory: Path = Path.cwd()
    vmax_intensity = np.percentile(matrix_data, percentile)
    plt.figure()
    color_map = LinearSegmentedColormap.from_list(
        "HiC", [(0, "white"), (0.5, "orange"), (1, "red")]
    )
    plt.imshow(matrix_data, cmap=color_map, vmin=0, vmax=vmax_intensity)
    # Create color bar for heatmap
    plt.colorbar()
    # Create dashed rectangles around genomic bin locations if present
    ax = plt.gca()
    tad1_rec = Rectangle(
        (tad1_genomic_data_start, tad1_genomic_data_start),
        tad1_length,
        tad1_length,
        fill=False,
        facecolor=None,
        linestyle="dashed",
    )
    tad2_rec = Rectangle(
        (tad2_genomic_data_start, tad2_genomic_data_start),
        tad2_length,
        tad2_length,
        fill=False,
        facecolor=None,
        linestyle="dashed",
    )
    ax.add_patch(tad1_rec)
    ax.add_patch(tad2_rec)
    plt.xticks(np.linspace(0, len(genomic_bin_matrix), 10), final_labels, rotation=45)
    plt.yticks(np.linspace(0, len(genomic_bin_matrix), 10), final_labels, rotation=45)
    plt.xlabel("Bin Region")
    plt.ylabel("Bin Region")
    plt.title(f"{title}")
    plt.savefig(f"{current_directory}/{filename}.png")
    plt.close()


## Part 2

# TODO: plot the GM matrix as a heatmap using the plot_hic_map() function
# TODO: plot the K562 matrix as a heatmap using the plot_hic_map() function

plot_hic_map(
    tad_genomic_bin_data_indices,
    gm_matrix,
    "GM HiC Heatmap",
    "gm_interaction_map",
    genomic_bin_matrix,
)
plot_hic_map(
    tad_genomic_bin_data_indices,
    k562_matrix,
    "K562 HiC Heatmap",
    "k562_interaction_map",
    genomic_bin_matrix,
)
## Part 3

# TODO: define a function block_contact_freq() that will subset the original matrix provided to it based on a list of index values provided to it.
# As input this function should accept an input matrix and index values which may be are the start and end locations of a TAD in that matrix
# The function should first generate a copy of the original matrix to work from (if using numpy we recommend np.copy() and if using pandas we recommend <variable_for_your_dataframe.copy()) <- will cover why in lecture
# The function should then slice the copied matrix at the positions indicated in the input list.
# Throw away the diagonal since this represents self-interactions within bins and is generally uninformative (set the diagoal values of the matrix to zero)
# calculate either the intra TAD interaction frequency if just a single TAD is in teh coordinates or the inter TAD interaction frequency if more than one TAD is present
# intra TAD interaction frequency formula: (sum_of_everything_in_subset_matrix)/(number_of_tads_in_original_matrix)
# inter TAD interaction frequency formula: (sum_of_everything_in_subset_matrix)/(number_of_tads_in_original_matrix) - (sum of each TAD's intra TAD interaction frequency)
# return the intra or inter TAD interaction frequncy


## Part 4

# TODO: use the block_contact_freq() function to calculate the intra-TAD interaction frequency for each TAD and the inter-TAD interaction frequency for GM12878 cells

# TODO: use the function to calculate the intra-TAD interaction frequency for each TAD and the inter-TAD interaction frequency for K562 cells

# TODO: convert all interaction frequencies into percents for each cell line. (if you are using pandas to store the interaction frequencies then make a copy of the dataframe you are using before doing so)

# TODO: Write all interaction frequencies intra and inter to the output file (cell_interaction_frequencies.txt). Then write each of the percents to the same outputfile.
