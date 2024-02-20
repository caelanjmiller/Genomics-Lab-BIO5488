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
        # Determine delimiter & if header exists for provided file
        tabular_delimiter = csv.Sniffer().sniff(rna_seq.read(5000), delimiters="\t").delimiter
        tabular_header: bool = csv.Sniffer().has_header(rna_seq.read(1024))
        file_contents: list = csv.reader(rna_seq, delimiter=tabular_delimiter)
        if tabular_header:
            next(rna_seq, None)
            return file_contents
        else:
            pass
        pass

def counts_per_million() -> dict:
    """"""
    pass


def fishers_linear_discriminant():
    """"""
    pass

RNA_SEQ_IO(EXPRESSION_DATA)

if len(argv) != 2:
    print(__doc__)
    exit(1)