#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
import matplotlib.pyplot as plt
import os
import csv

"""
Python script to parse CGI average methylation .bed file & return:
- Plot of average methylation across all CGIs
Usage: python3 analyze_CGI_methylation.py <CGI AVERAGE METHYLATION BED FILE>
"""

CGI_METHYLATION_BED = Path(argv[1])


class CGI:
    def __init__(
        self,
        chromosome_number: str,
        start_coordinate: int,
        stop_coordinate: int,
        cpg_name: str,
        average_methylation_level: float,
    ):
        self.chromosome_number = chromosome_number
        self.start_coordinate = start_coordinate
        self.stop_coordinate = stop_coordinate
        self.cpg_name = cpg_name
        self.average_methylation_level = average_methylation_level

    def __repr__(self):
        return f"CGI({self.chromosome_number}, {self.start_coordinate}, {self.stop_coordinate}, {self.cpg_name}, {self.average_methylation_level})"


def BED_IO(FILE) -> list:
    """Parse in BED file & instantiate a CGI object for each entry (row)"""
    extensions: str = ".bed"
    bed_filename: str = os.path.basename(FILE)
    if bed_filename.endswith(extensions):
        try:
            with open(FILE) as bed:
                file_contents: list = csv.reader(bed, delimiter="\t")
                cgis: list = []
                for element in file_contents:
                    chromosome_number: str = str(element[0])
                    start_coordinate: int = int(element[1])
                    stop_coordinate: int = int(element[2])
                    cpg_name: str = element[3]
                    average_methylation_level: float = float(element[4])
                    cgi: CGI = CGI(
                        chromosome_number,
                        start_coordinate,
                        stop_coordinate,
                        cpg_name,
                        average_methylation_level,
                    )
                    cgis.append(cgi)
            return cgis
        except Exception as e:
            return f"{e}"
    else:
        return f"Unable to open {bed_filename}"


def create_CpG_methylation_distribution(cgis: list, FILE):
    """Create histogram of CGI methylation levels"""
    current_directory: Path = Path.cwd()
    basename: str = os.path.basename(FILE).split(".bed")[0]
    CGI_methylation: list = [cgi.average_methylation_level for cgi in cgis]
    plt.hist(CGI_methylation, bins=20)
    plt.xlabel("Average CGI Methylation Levels")
    plt.ylabel("Frequency")
    plt.title(f"CGI Methylation Levels in {cgis[0].chromosome_number}")
    plt.savefig(f"{current_directory}/{basename}_distribution.png")


cgis: list = BED_IO(CGI_METHYLATION_BED)
create_CpG_methylation_distribution(cgis, CGI_METHYLATION_BED)

if len(argv) != 2:
    print(__doc__)
    exit(1)
