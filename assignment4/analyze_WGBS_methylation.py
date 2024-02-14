#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
import matplotlib.pyplot as plt
import os
import csv

WBGS_BED = Path(argv[1])
PRINTOUT = argv[2]


# Create a class for each row in BED file
class Bed_Row:
    def __init__(
        self,
        chromosome_number: int,
        start_coordinate: int,
        stop_coordinate: int,
        as_c: int,
        as_t: int,
    ):
        self.chromosome_number = chromosome_number
        self.start_coordinate = start_coordinate
        self.stop_coordinate = stop_coordinate
        self.as_c = as_c
        self.as_t = as_t

    def calculate_CpG_methylation_levels(self):
        """Calculate CpG methylation levels from cytosine & thymine ratios"""
        self.methylation_level: float = round((self.as_c / (self.as_c + self.as_t)), 3)

    def __str__(self) -> str:
        return (
            f"Chromosome: {self.chromosome_number}\n"
            f"Start Coordinate: {self.start_coordinate}\n"
            f"Stop Coordinate: {self.stop_coordinate}\n"
            f"Times base called as C: {self.as_c}\n"
            f"Times base called as T: {self.as_t}\n"
        )


def BED_IO(FILE) -> list:
    extensions: str = ".bed"
    bed_filename: str = os.path.basename(FILE)
    if bed_filename.endswith(extensions):
        try:
            with open(FILE) as bed:
                file_contents: list = csv.reader(bed, delimiter="\t")
                bed_coordinates: list = []
                for element in file_contents:
                    bed_row: Bed_Row = Bed_Row(
                        element[0],
                        element[1],
                        element[2],
                        element[3],
                        element[4],
                    )
                    bed_coordinates.append(bed_row)
            return bed_coordinates
        except Exception:
            return f"Unable to open {bed_filename}"


def assign_methylation_levels(bed_file: list):
    """Calculate methylation levels for each of the rows of a WGBS BED file"""
    for row in bed_file:
        row.calculate_CpG_methylation_levels()


def output_bed_file(bed_file: list, FILE):
    """Output BED file (tab-delimited) with CpG methylation - excludes 0X coverage"""
    # Get current directory
    current_directory: Path = Path.cwd()
    # Extract basename of BED file to use for new methylation bed file
    basename: str = os.path.basename(FILE).split(".bed")[0]
    with open(f"{current_directory}/{basename}_CpG_methylation.bed") as bed:
        # Exclude rows with 0X coverage
        for row in bed_file:
            if row.methylation == 0:
                continue
            else:
                bed.write(
                    f"{row.chromosome_number}\t{row.start_coordinate}\t{row.stop_coordinate}\t{row.methylation}\n"
                )


if len(argv) != 2:
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

bed_file: list = BED_IO(WBGS_BED)

if PRINTOUT == "methylation-bed-file":
    output_bed_file(bed_file)
elif PRINTOUT == "methylation-level-plot":
    pass
elif PRINTOUT == "read-coverage-plot":
    pass
elif PRINTOUT == "zero-coverage":
    pass
else:
    raise Exception(
        f"Provide valid printout option:\n"
        f"methylation-bed-file : Create bed file of CpG methylation levels\n"
        f"methylation-level-plot : Create histogram of CpG methylation levels"
        f"read-coverage-plot : Create histogram of CpG read coverage\n"
        f"zero-coverage : Prints out fraction of CpGs with 0X read coverage\n"
    )
