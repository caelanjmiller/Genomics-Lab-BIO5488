#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import csv
import os
""""
Python script to 
"""

EXPRESSION_DATA = Path(argv[1])


def RNA_SEQ_IO(FILE: Path) -> dict:
    """
    Parse in RNA Seq Data and return a dict with gene
    as the key and sample counts as the value
    """
    with open(FILE, "r") as rna_seq:
        # Determine if header exists for provided file
        file_contents: list = csv.reader(rna_seq, delimiter="\t")
        rna_seq_data: dict = {}
        next(rna_seq, None)
        for gene_count_data in file_contents:
            gene_name: str = gene_count_data[0]
            # Convert counts from str to int
            rna_seq_data[gene_name] = [*map(int, gene_count_data[1:])]
        # Return sorted dict
        return dict(sorted(rna_seq_data.items()))


def add_sample_number(rna_seq_data: dict) -> dict:
    """Add Sample Name (e.g. Before_1 ) to list of gene count data"""
    # Initialize defaultdict to create nested dictionary - [Gene][Sample_Name] = [Count_Data]
    annotated_rna_seq_data: defaultdict = defaultdict(dict)
    number_samples: int = int((len(next(iter(rna_seq_data.values()))) / 2))
    for gene_name, count_data in rna_seq_data.items():
        for index, sample_count in enumerate(count_data):
            if index <= number_samples - 1:
                sample_name: str = f"Before_{index + 1}"
                annotated_rna_seq_data[gene_name][sample_name] = sample_count
            else:
                sample_name: str = f"After_{(index - number_samples) + 1}"
                annotated_rna_seq_data[gene_name][sample_name] = sample_count
    return dict(sorted(annotated_rna_seq_data.items()))


def transpose_data(rna_seq_data: dict) -> dict:
    """
    Transpose data from:
    Gene: Gene Counts - Sample
    Sample: Counts Per Gene
    """
    # Initialize a default dict to return a dict with structure: [Sample][Gene] = [Count]
    # e.g. {'Before_1': {'A2ML1': 4}}
    transposed_rna_seq_data: defaultdict = defaultdict(dict)
    for gene_name, count_data in rna_seq_data.items():
        for sample_name, count in count_data.items():
            transposed_rna_seq_data[sample_name][gene_name] = count
    return dict(transposed_rna_seq_data.items())


def filter_zero_gene_expression(rna_seq_data: dict) -> dict:
    """Filter out genes with 0 count data within a RNA Seq dataset"""
    filtered_rna_seq_data: dict = {}
    for gene, count_data in rna_seq_data.items():
            # If sum of list is greater than 0, append to new dict
            if sum(count_data) > 0:
                filtered_rna_seq_data[gene] = count_data
            else:
                continue
    return dict(sorted(filtered_rna_seq_data.items()))


def calculate_library_sizes(transposed_rna_seq_data: dict) -> dict:
    """Calculate library size of each RNA Seq sample"""
    library_sizes: dict = {}
    for sample_name, count_data in transposed_rna_seq_data.items():
        library_sizes[sample_name] = sum(count_data.values())
    return library_sizes


def counts_per_million(filtered_rna_seq_data: dict, library_sizes: dict) -> dict:
    """Calculate counts per million for each gene"""
    rna_seq_cpm: dict = {}
    total_library_sizes: list = list(library_sizes.values())
    for gene_name, count_data in filtered_rna_seq_data.items():
            rna_seq_cpm[gene_name] = [((10**6) * (raw_count / total_counts))for raw_count, total_counts in zip(count_data, total_library_sizes)]
    return dict(sorted(rna_seq_cpm.items()))
    
def filter_by_counts_per_million(filtered_rna_seq_data: dict, rna_seq_cpm: dict, threshold: float, num_fail: int) -> dict:
    """Filter RNA Seq data dict by CPM threshold"""
    cpm_filtered_rna_seq_data: dict = {}
    for gene_name, cpm_data in rna_seq_cpm.items():
        number_samples_fail_cpm: list = []
        for cpm in cpm_data:
            if cpm < threshold:
                number_samples_fail_cpm.append(cpm)
        if len(number_samples_fail_cpm) < num_fail:
           cpm_filtered_rna_seq_data[gene_name] = filtered_rna_seq_data[gene_name]
    return dict(sorted(cpm_filtered_rna_seq_data.items()))

def calculate_library_size_range(rna_seq_data: dict) -> tuple:
    """Calculate range of library size from RNA Seq sample data"""
    pass


def create_library_size_histogram(cpm_filtered_rna_seq: dict) -> None:
    """Create histogram of library sizes for each sample in RNA Seq dataset"""
    current_directory: Path = Path.cwd()
    annotated_rna_seq_data: dict = add_sample_number(cpm_filtered_rna_seq)
    transposed_cpm_filtered_rna_seq_data: dict = transpose_data(annotated_rna_seq_data)
    total_sample_library_sizes: dict = calculate_library_sizes(transposed_cpm_filtered_rna_seq_data)
    plt.bar(x=total_sample_library_sizes.keys(), height=total_sample_library_sizes.values())
    plt.xlabel("Samples")
    plt.ylabel("Library Size (in Millions)")
    plt.title(f"Library Sizes of RNA Seq Samples")
    plt.savefig(f"{current_directory}/library_size.png")

def fishers_linear_discriminant():
    """"""
    pass


rna_seq_data: dict = RNA_SEQ_IO(EXPRESSION_DATA)
# print(f"Number of Genes: {len(rna_seq_data.keys())}")
filtered_rna_seq_data: dict = filter_zero_gene_expression(rna_seq_data)
# print(f"Number of Genes After Filtering: {len(filtered_rna_seq_data.keys())}")
annotated_rna_seq_data: dict = add_sample_number(filtered_rna_seq_data)
transposed_unfiltered_rna_seq_data: dict = transpose_data(
    annotated_rna_seq_data
)
rna_seq_library_sizes: dict = calculate_library_sizes(transposed_unfiltered_rna_seq_data)
rna_seq_cpm: dict = counts_per_million(filtered_rna_seq_data, rna_seq_library_sizes)
cpm_filtered_rna_seq: dict = filter_by_counts_per_million(filtered_rna_seq_data,rna_seq_cpm, float(1), 20)
# print(f"Number of genes left after CPM filtration: {len(cpm_filtered_rna_seq.keys())}")
create_library_size_histogram(cpm_filtered_rna_seq)
if len(argv) != 2:
    print(__doc__)
    exit(1)
