from pathlib import Path
from collections import Counter
from itertools import zip_longest
from sys import argv, exit
import csv

"""
Python script to demultiplex scRNA Seq cell barcodes
Usage: python3 demultiplex.py <CALL MATRIX CSV>
"""

if len(argv) != 2:
    print(__doc__)
    exit(1)


CALL_MATRIX = Path(argv[1])


class Cell:
    def __init__(
        self,
        barcode: str,
        control_count: tuple,
        treatment_count: tuple,
        cell_tag1_count: tuple,
        cell_tag2_count: tuple,
    ) -> None:
        self.barcode = barcode
        self.control_count = control_count
        self.treatment_count = treatment_count
        self.cell_tag1_count = cell_tag1_count
        self.cell_tag2_count = cell_tag2_count
        self.set_cell_identity()

    def set_cell_identity(self) -> None:
        """Set cellular identity based upon bar code count"""
        if self.control_count[1] == 1:
            self.identity = "control"
        elif self.treatment_count[1] == 1:
            self.identity = "treatment"
        else:
            self.identity = "undetermined"


def parse_cell_matrix_csv(FILE: Path) -> list:
    """Parse cell matrix CSV & initialize Cell objects from provided information"""
    with open(FILE, "r") as CSV:
        scrna_seq_data: list = []
        file_contents: list = csv.reader(CSV, delimiter=",")
        csv_header: list = next(file_contents)
        for line in file_contents:
            cell: Cell = Cell(
                line[0],
                tuple((csv_header[1], int(line[1]))),
                tuple((csv_header[2], int(line[2]))),
                tuple((csv_header[3], int(line[3]))),
                tuple((csv_header[4], int(line[4]))),
            )
            scrna_seq_data.append(cell)
    return scrna_seq_data


def initial_demultiplexing(scrna_seq_data: list) -> dict:
    """Parse scRNA Seq cell barcoding to assign cells to treatment, control or nondetermined groups"""
    cell_identity_counts: dict = dict(
        Counter([cell.identity for cell in scrna_seq_data])
    )
    return cell_identity_counts


def secondary_demultiplexing(scrna_seq_data: list) -> dict:
    """
    Parse scRNA Seq cell barcoding to assign cells to treatment, control or nondetermined groups
    after correcting for a Hamming Distance of 1
    """
    pass


def calculate_hamming_distance(sequence_one: str, sequence_two: str) -> int:
    """Given two strings, return the Hamming Distance between them"""
    hamming_distance: int = 0
    for i, j in zip_longest(list(sequence_one), list(sequence_two)):
        if i != j:
            hamming_distance += 1
        else:
            hamming_distance += 0
    return hamming_distance


scrna_seq_data: list = parse_cell_matrix_csv(CALL_MATRIX)
cell_identity_counts: dict = initial_demultiplexing(scrna_seq_data)

print("....Initial Demultiplexing....")
for cell_identity, count in cell_identity_counts.items():
    print(f"{cell_identity}: {count}")
