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
        chromosome_number: str,
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
        if self.as_c == 0:
            self.methylation_level: float = float(0)
        else:
            self.methylation_level: float = round(
                (self.as_c / (self.as_c + self.as_t)), 3
            )

    def calculate_read_coverage(self):
        """Calculate coverage of reads for CpGs at a given genomic location"""
        self.read_coverage: int = self.as_c + self.as_t

    def __repr__(self):
        return f"Bed_Row({self.chromosome_number}, {self.start_coordinate}, {self.stop_coordinate}, {self.as_c}, {self.as_t}, {self.methylation_level}, {self.base_coverage})"

    def __str__(self):
        return (
            f"Chromosome: {self.chromosome_number}\n"
            f"Start Coordinate: {self.start_coordinate}\n"
            f"Stop Coordinate: {self.stop_coordinate}\n"
            f"Times base called as C: {self.as_c}\n"
            f"Times base called as T: {self.as_t}\n"
            f"Read Coverage: {self.base_coverage}\n"
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
                    chromosome_number: str = str(element[0])
                    start_coordinate: int = int(element[1])
                    stop_coordinate: int = int(element[2])
                    as_c: int = int(element[3])
                    as_t: int = int(element[4])
                    bed_row: Bed_Row = Bed_Row(
                        chromosome_number, start_coordinate, stop_coordinate, as_c, as_t
                    )
                    bed_coordinates.append(bed_row)
            return bed_coordinates
        except Exception as e:
            return f"{e}"
    else:
        return f"Unable to open {bed_filename}"


def assign_methylation_levels(bed_coordinates: list):
    """Assign methylation levels for each of the rows of a WGBS BED file"""
    for bed_row in bed_coordinates:
        bed_row.calculate_CpG_methylation_levels()


def assign_read_coverage(bed_coordinates: list):
    """Assign read coverage for each of the rows of WGBS BED file"""
    for bed_row in bed_coordinates:
        bed_row.calculate_read_coverage()


def output_CpG_bed_file(bed_coordinates: list, FILE):
    """Output BED file (tab-delimited) with CpG methylation - excludes 0X coverage"""
    # Get current directory
    current_directory: Path = Path.cwd()
    # Extract basename of BED file to use for new methylation bed file
    basename: str = os.path.basename(FILE).split(".bed")[0]
    with open(f"{current_directory}/{basename}_CpG_methylation.bed", "w") as bed:
        # Exclude rows with 0X coverage
        for bed_row in bed_coordinates:
            if bed_row.methylation_level == 0:
                continue
            else:
                bed.write(
                    f"{bed_row.chromosome_number}\t{bed_row.start_coordinate}\t{bed_row.stop_coordinate}\t{bed_row.methylation_level}\n"
                )


def create_CpG_methylation_distribution(bed_coordinates: list, FILE):
    """Create histogram of CpG methylation levels"""
    current_directory: Path = Path.cwd()
    basename: str = os.path.basename(FILE).split(".bed")[0]
    CpG_methylation_nonzero: list = [
        bed_row.methylation_level for bed_row in bed_coordinates
    ]
    plt.hist(CpG_methylation_nonzero, bins=20)
    plt.xlabel("CpG Methylation Levels")
    plt.ylabel("Frequency")
    plt.title(f"CpG Methylation Levels in {basename}")
    plt.savefig(f"{current_directory}/{basename}_methylation_distribution.png")


def create_CpG_read_coverage_distribution(bed_coordinates: list, FILE):
    """Create histogram of CpG read coverage"""
    current_directory: Path = Path.cwd()
    basename: str = os.path.basename(FILE).split(".bed")[0]
    CpG_read_coverage: list = [
        bed_row.read_coverage
        for bed_row in bed_coordinates
        if bed_row.read_coverage in range(0, 101)
    ]
    plt.hist(CpG_read_coverage, bins=20)
    plt.xlabel("CpG Read Coverage")
    plt.ylabel("Frequency")
    plt.title(f"CpG Read Coverage in {basename}")
    plt.savefig(f"{current_directory}/{basename}_CpG_coverage_distribution.png")


bed_coordinates: list = BED_IO(WBGS_BED)
assign_methylation_levels(bed_coordinates)
assign_read_coverage(bed_coordinates)

if PRINTOUT == "methylation-bed-file":
    output_CpG_bed_file(bed_coordinates, WBGS_BED)
elif PRINTOUT == "methylation-level-plot":
    create_CpG_methylation_distribution(bed_coordinates, WBGS_BED)
elif PRINTOUT == "read-coverage-plot":
    create_CpG_read_coverage_distribution(bed_coordinates, WBGS_BED)
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

if len(argv) != 3:
    docstring = (
        f"Python script to parse methylation .bed files & returns:\n"
        f"1. BED File of CpG methylation for a user provided whole genome bisulfite sequence (WGBS) BED file\n"
        f"2. Plots for distribution of CpG methylation levels\n"
        f"3. Plots for distribution of read coverage for all CpGs [0-100X]\n"
        f"4. Print CpG fraction with 0X read coverage\n"
        f"Usage: python3 analyze_WGBS_methylation.py <WGBS BED FILE> <PRINTOUT>"
    )
    print(docstring)
    exit(1)
