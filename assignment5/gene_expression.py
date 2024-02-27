#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
import matplotlib.pyplot as plt
from statistics import mean
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

def transpose_data(rna_seq_data: dict) -> dict:
    """
    Transpose data from:
    Gene: Gene Counts - Sample
    Sample: Counts Per Gene
    """
    gene_names: list = []
    for gene, count_data in rna_seq_data.items():
        number_samples: int = (len(count_data) / 2)

    pass

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

if len(argv) != 2:
    print(__doc__)
    exit(1)