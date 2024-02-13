#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
import matplotlib.pyplot as plt
from ast import literal_eval
import os
import csv

WBGS_BED = Path(argv[1])
CHR_HEADER = argv[2]
PRINTOUT = argv[3]


def BED_IO(FILE):
    extensions: str = ".bed"
    bed_filename: str = os.path.basename(FILE)
    if bed_filename.endswith(extensions):
        try:
            with open(FILE) as bed:
                file_contents: list = csv.reader(bed, delimiter="\t")
                bed_as_dict: dict = {}
                for element in file_contents:
                    file_string = element[0].split("\t")
                    file_key: str = file_string[0]
                    # Use literal_eval to interpret string as a dictionary
                    file_data: dict = literal_eval(file_string[1])
                    bed_as_dict[file_key] = file_data
        except Exception:
            return f"Unable to open {bed_filename}"
        return bed_as_dict


def calculate_CpG_methylation_levels(bed: dict) -> dict:
    """Calculate CpG methylation levels from a whole genome bisulfite sequencing BED file"""
    pass


if len(argv) != 4:
    docstring = (
        f"Python script to parse methylation .bed files & returns:\n"
        f"1. BED File of CpG methylation for a user provided whole genome bisulfite sequence (WGBS) BED file\n"
        f"2. Plots for distribution of CpG methylation levels\n"
        f"3. Plots for distribution of read coverage for all CpGs [0-100X]\n"
        f"4. Print CpG fraction with 0X read coverage\n"
        f"Usage: python3 analyze_WGBS_methylation.py <WGBS BED FILE> <CHR HEADER> <PRINTOUT>"
    )
    print(docstring)
    exit(1)
