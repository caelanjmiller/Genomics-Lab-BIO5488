#!/usr/bin/env python3

from sys import argv, exit
from pathlib import Path
import os
import csv

"""
Python script that takes BED file of gene coordinates and returns
BED file of their promoters
Usage: python3 generate_promoters.py <BED OF COORDINATES>
"""

GENE_COORDINATES_BED = Path(argv[1])


class Gene:
    def __init__(
        self,
        chromosome_number: str,
        start_coordinate: int,
        stop_coordinate: int,
        gene_name: str,
        number_exons: int,
        dna_strand: str,
    ):
        self.chromosome_number = chromosome_number
        self.start_coordinate = start_coordinate
        self.stop_coordinate = stop_coordinate
        self.gene_name = gene_name
        self.number_exons = number_exons
        self.dna_strand = dna_strand

    def generate_promoter_regions(self):
        """Generate promoter regions utilizing gene start & stop coordinates [1st 1000bp as promoter]"""
        if self.dna_strand == "+":
            self.promoter_region: tuple = tuple(
                (self.start_coordinate, self.start_coordinate + 1000)
            )
        elif self.dna_strand == "-":
            self.promoter_region: tuple = tuple(
                (self.stop_coordinate - 1000, self.stop_coordinate)
            )

    def __repr__(self):
        return f"CGI({self.chromosome_number}, {self.start_coordinate}, {self.stop_coordinate}, {self.gene_name}, {self.promoter_region}, {self.number_exons}, {self.dna_strand})"


def BED_IO(FILE) -> list:
    """Parse in BED file & instantiate a CGI object for each entry (row)"""
    extensions: str = ".bed"
    bed_filename: str = os.path.basename(FILE)
    if bed_filename.endswith(extensions):
        try:
            with open(FILE) as bed:
                file_contents: list = csv.reader(bed, delimiter="\t")
                genes: list = []
                for element in file_contents:
                    chromosome_number: str = str(element[0])
                    start_coordinate: int = int(element[1])
                    stop_coordinate: int = int(element[2])
                    gene_name: str = element[3]
                    number_exons: int = int(element[4])
                    dna_strand: str = element[5]
                    gene: Gene = Gene(
                        chromosome_number,
                        start_coordinate,
                        stop_coordinate,
                        gene_name,
                        number_exons,
                        dna_strand,
                    )
                    genes.append(gene)
            return genes
        except Exception as e:
            return f"{e}"
    else:
        return f"Unable to open {bed_filename}"


def assign_promoter_regions(genes: list):
    """Assign promoter regions for each Gene object"""
    for gene in genes:
        gene.generate_promoter_regions()


def output_gene_promoters_bed_file(genes: list, FILE):
    """Output BED file (tab-delimited) with gene information"""
    # Get current directory
    current_directory: Path = Path.cwd()
    # Extract basename of BED file to use for new bed file
    basename: str = os.path.basename(FILE).split(".bed")[0]
    with open(f"{current_directory}/{basename}_promoters.bed", "w") as bed:
        # Write information
        for gene in genes:
            promoter_start_coordinate, promoter_stop_coordinate = gene.promoter_region
            bed.write(
                f"{gene.chromosome_number}\t{promoter_start_coordinate}\t{promoter_stop_coordinate}\t{gene.gene_name}\t{gene.dna_strand}\n"
            )


genes: list = BED_IO(GENE_COORDINATES_BED)
assign_promoter_regions(genes)
output_gene_promoters_bed_file(genes, GENE_COORDINATES_BED)

if len(argv) != 2:
    print(__doc__)
    exit(1)
