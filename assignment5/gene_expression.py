#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
from collections import defaultdict
import csv

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
        tabular_header: bool = csv.Sniffer().has_header(rna_seq.read(1024))
        file_contents: list = csv.reader(rna_seq, delimiter="\t")
        if tabular_header:
            rna_seq_data: dict = {}
            next(rna_seq, None)
            for gene_count_data in file_contents:
                gene_name: str = gene_count_data[0]
                # Convert counts from str to int
                rna_seq_data[gene_name] = [*map(int, gene_count_data[1:])]
            # Return sorted dict
            return dict(sorted(rna_seq_data.items()))
        else:
            rna_seq_data: dict = {}
            for gene_count_data in file_contents:
                gene_name: str = gene_count_data[0]
                rna_seq_data[gene_name] = [*map(int, gene_count_data[1:])]
            return dict(sorted(rna_seq_data.items()))


def create_transposed_sample_dict(rna_seq_data: dict) -> dict:
    """Create skeleton dict for use in transposing data"""
    # Grab an arbitrary dict element and calculate the number of unique Before & After Sample number
    number_samples: int = int((len(next(iter(rna_seq_data.values()))) / 2))
    transposed_skeleton_data: dict = {}
    for sample_number in range(1, number_samples + 1):
        transposed_skeleton_data[f"Before_{sample_number}"] = {}
    for sample_number in range(1, number_samples + 1):
        transposed_skeleton_data[f"After_{sample_number}"] = {}
    return transposed_skeleton_data


def add_sample_number(rna_seq_data: dict) -> dict:
    """Add Sample Name (e.g. Before_1 ) to list of gene count data"""
    # Initialize defaultdict to create nested dictionary - [Gene][Sample_Name] = [Count_Data]
    annotated_rna_seq_data: defaultdict = defaultdict(dict)
    number_samples: int = int((len(next(iter(rna_seq_data.values()))) / 2))
    for gene_name, count_data in rna_seq_data.items():
        for index, sample_count in enumerate(count_data):
            if index <= number_samples:
                sample_name: str = f"Before_{index}"
                annotated_rna_seq_data[gene_name][sample_name] = sample_count
            else:
                sample_name: str = f"After_{index - number_samples}"
                annotated_rna_seq_data[gene_name][sample_name] = sample_count
    return dict(sorted(annotated_rna_seq_data.items()))


def transpose_data(rna_seq_data: dict, transposed_skeleton_data: dict) -> dict:
    """
    Transpose data from:
    Gene: Gene Counts - Sample
    Sample: Counts Per Gene
    """
    sample_names: list = transposed_skeleton_data.keys()
    for gene_name, count_data in rna_seq_data.items():
        for count in count_data:
            pass
    return sample_names


def filter_zero_gene_expression(transposed_data: dict) -> dict:
    """Filter out genes with 0 count data within a RNA Seq dataset"""
    pass


def counts_per_million() -> dict:
    """"""
    pass


def fishers_linear_discriminant():
    """"""
    pass


rna_seq_data: dict = RNA_SEQ_IO(EXPRESSION_DATA)
annotated_rna_seq_data: dict = add_sample_number(rna_seq_data)
transposed_skeleton_data: dict = create_transposed_sample_dict(rna_seq_data)
# transposed_rna_seq_data = transpose_data(rna_seq_data, transposed_skeleton_data)

if len(argv) != 2:
    print(__doc__)
    exit(1)
