#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
from collections import defaultdict
from itertools import islice
import matplotlib.pyplot as plt
import numpy as np
import csv

""""
Python script to analyze RNA Seq Data
Usage: python3 gene_expression.py <RNA SEQ DATA>
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
        rna_seq_cpm[gene_name] = [
            ((10**6) * (raw_count / total_counts))
            for raw_count, total_counts in zip(count_data, total_library_sizes)
        ]
    return dict(sorted(rna_seq_cpm.items()))


def filter_by_counts_per_million(
    filtered_rna_seq_data: dict, rna_seq_cpm: dict, threshold: float, num_fail: int
) -> dict:
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
    library_sizes: list = list(rna_seq_data.values())
    # Round to nearest integer (cognizant of rounding to the closest int)
    range: tuple = tuple((round(min(library_sizes)), round((max(library_sizes)))))
    return range


def create_library_size_histogram(cpm_filtered_rna_seq: dict) -> None:
    """Create histogram of library sizes for each sample in RNA Seq dataset"""
    current_directory: Path = Path.cwd()
    annotated_rna_seq_data: dict = add_sample_number(cpm_filtered_rna_seq)
    transposed_cpm_filtered_rna_seq_data: dict = transpose_data(annotated_rna_seq_data)
    total_sample_library_sizes: dict = calculate_library_sizes(
        transposed_cpm_filtered_rna_seq_data
    )
    # Scale library sizes for easier visualization
    scaled_library_sizes: list = list(
        map(
            lambda library_size: library_size / 10**6,
            list(total_sample_library_sizes.values()),
        )
    )
    plt.figure(figsize=(8, 6))
    plt.bar(
        x=total_sample_library_sizes.keys(),
        height=(scaled_library_sizes),
        color=np.random.rand(len(scaled_library_sizes), 3),
    )
    plt.xlabel("Samples")
    plt.ylabel("Library Size (in Millions)")
    plt.title("Library Sizes Prior to Normalization")
    plt.xticks(rotation=45, ha="right")
    plt.savefig(f"{current_directory}/library_size.png")
    plt.close()


def calculate_upper_quartile_normalization(cpm_filtered_rna_seq_data: dict) -> dict:
    """Calculate the upper quartile normalized count for RNA Seq samples"""
    upper_quartile_normalization_data: dict = {}
    # Create transposed dictionary
    transposed_annotated_cpm_filtered_rna_seq_data: dict = transpose_data(
        add_sample_number(cpm_filtered_rna_seq_data)
    )
    # Create list of upper quartile values for each sample
    upper_quartile_values: list = []
    for count_data in transposed_annotated_cpm_filtered_rna_seq_data.values():
        upper_quartile_values.append(np.percentile(list(count_data.values()), 75))
    # Calculate mean value of upper quartiles across all samples
    mean_value: float = np.mean(upper_quartile_values)
    for gene_name, count_data in cpm_filtered_rna_seq_data.items():
        # Perform upper quartile normalization
        upper_quartile_normalization_data[gene_name] = [
            ((raw_count / upper_quartile) * mean_value)
            for raw_count, upper_quartile in zip(count_data, upper_quartile_values)
        ]
    return upper_quartile_normalization_data


def create_library_size_normalized_histogram(
    transposed_annotated_normalized_rna_seq_data: dict,
) -> None:
    """Create histogram of normalized library sizes for each sample in RNA Seq dataset"""
    current_directory: Path = Path.cwd()
    total_sample_library_sizes: dict = calculate_library_sizes(
        transposed_annotated_normalized_rna_seq_data
    )
    # Scale library sizes for easier visualization
    scaled_library_sizes: list = list(
        map(
            lambda library_size: library_size / 10**6,
            list(total_sample_library_sizes.values()),
        )
    )
    plt.figure(figsize=(8, 6))
    plt.bar(
        x=total_sample_library_sizes.keys(),
        height=(scaled_library_sizes),
        color=np.random.rand(len(scaled_library_sizes), 3),
    )
    plt.xlabel("Samples")
    plt.ylabel("Library Size (in Millions)")
    plt.title("Library Sizes After Normalization")
    plt.xticks(rotation=45, ha="right")
    plt.savefig(f"{current_directory}/library_size_normalized.png")
    plt.close()


def calculate_fishers_linear_discriminant(
    filtered_normalized_rna_seq_data: dict,
) -> dict:
    """Calculate Fisher's Linear Discriminant (FLD) for determination of differential expression"""
    fishers_linear_discriminant_data: dict = {}
    number_samples: int = int(
        (len(next(iter(filtered_normalized_rna_seq_data.values()))) / 2)
    )
    for gene, count_data in filtered_normalized_rna_seq_data.items():
        before_count_data: list = count_data[
            number_samples - number_samples : number_samples
        ]
        after_count_data: list = count_data[number_samples:]
        fishers_linear_discriminant_data[gene] = (
            np.mean(before_count_data) - np.mean(after_count_data)
        ) ** 2 / ((np.std(before_count_data)) ** 2 + (np.std(after_count_data) ** 2))
    return fishers_linear_discriminant_data


def capture_top_genes_FLD(fishers_linear_discriminant_data: dict) -> dict:
    """Return top 10 genes with highest FLD score"""
    # Sort by FLD & grab top 10 highest FLD scores
    ascending_FLD_scores: dict = dict(
        (
            sorted(
                fishers_linear_discriminant_data.items(),
                key=lambda FLD_score: FLD_score[1],
                reverse=True,
            )
        )
    )
    top_FLD_scores: dict = dict(islice(ascending_FLD_scores.items(), 10))
    return top_FLD_scores


def print_top_FLD_scores(top_FLD_scores: dict) -> None:
    """Print out top FLD scores to stdout"""
    for gene, FLD_score in top_FLD_scores.items():
        print(f"{gene}: {FLD_score}")


def create_expression_level_bar_chart(
    upper_normalization_data: dict, gene: str
) -> None:
    """Create a bar graph of the mean expression of a supplied gene"""
    merged_data: dict = {}
    before_count_sem: float = 0.0
    after_count_sem: float = 0.0
    for gene_name, count_data in upper_normalization_data.items():
        number_samples: int = int(
            (len(next(iter(upper_normalization_data.values()))) / 2)
        )
        if gene == gene_name:
            # Assign Before & After Gene Sample normalized count data
            before_count_data: list = count_data[
                number_samples - number_samples : number_samples
            ]
            after_count_data: list = count_data[number_samples:]
            merged_data["Before"] = np.mean(before_count_data)
            merged_data["After"] = np.mean(after_count_data)
            # Calculation of standard error of mean (SEM)
            before_count_sem = np.std(before_count_data) / np.sqrt(
                np.size(before_count_data)
            )
            after_count_sem = np.std(after_count_data) / np.sqrt(
                np.size(after_count_data)
            )
    current_directory: Path = Path.cwd()
    plt.figure(figsize=(8, 6))
    plt.bar(
        ["Before", "After"],
        list(merged_data.values()),
        color=np.random.rand(2, 3),
    )
    plt.errorbar(
        ["Before", "After"],
        list(merged_data.values()),
        yerr=[before_count_sem, after_count_sem],
        fmt="o",
        color="#000000",
    )
    plt.xlabel("Group")
    plt.ylabel("Mean Expression Level")
    plt.title(f"Mean Expression Levels of {gene.upper()}")
    plt.xticks(ha="right")
    plt.savefig(f"{current_directory}/mean_expression.png")
    plt.close()


rna_seq_data: dict = RNA_SEQ_IO(EXPRESSION_DATA)
print(f"Number of Genes: {len(rna_seq_data.keys())}")
zero_filtered_rna_seq_data: dict = filter_zero_gene_expression(rna_seq_data)
print(
    f"Number of genes after filtering zero expression: {len(zero_filtered_rna_seq_data.keys())}"
)
annotated_rna_seq_data: dict = add_sample_number(zero_filtered_rna_seq_data)
transposed_unfiltered_rna_seq_data: dict = transpose_data(annotated_rna_seq_data)
rna_seq_library_sizes: dict = calculate_library_sizes(
    transposed_unfiltered_rna_seq_data
)
rna_seq_cpm: dict = counts_per_million(
    zero_filtered_rna_seq_data, rna_seq_library_sizes
)
cpm_filtered_rna_seq_data: dict = filter_by_counts_per_million(
    zero_filtered_rna_seq_data, rna_seq_cpm, float(1), 20
)
transposed_annotated_cpm_filtered_rna_seq_data: dict = transpose_data(
    add_sample_number(cpm_filtered_rna_seq_data)
)
print(
    f"Number of genes left after CPM filtration: {len(cpm_filtered_rna_seq_data.keys())}"
)
print(
    f"Range of Library Sizes: {calculate_library_size_range(calculate_library_sizes(transposed_annotated_cpm_filtered_rna_seq_data))}"
)
create_library_size_histogram(cpm_filtered_rna_seq_data)
upper_normalization_data: dict = calculate_upper_quartile_normalization(
    cpm_filtered_rna_seq_data
)
transposed_annotated_normalized_rna_seq_data: dict = transpose_data(
    add_sample_number(upper_normalization_data)
)
create_library_size_normalized_histogram(transposed_annotated_normalized_rna_seq_data)
print(
    f"Range of Library Sizes: {calculate_library_size_range(calculate_library_sizes(transposed_annotated_normalized_rna_seq_data))}"
)
fishers_linear_discriminant_data: dict = calculate_fishers_linear_discriminant(
    upper_normalization_data
)
print(fishers_linear_discriminant_data)
print_top_FLD_scores(capture_top_genes_FLD(fishers_linear_discriminant_data))
create_expression_level_bar_chart(upper_normalization_data, "FNDC5")
if len(argv) != 2:
    print(__doc__)
    exit(1)
# Caelan Miller - 2024
