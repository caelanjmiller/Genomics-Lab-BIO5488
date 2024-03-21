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

    def update_cell_identity(self, identity: str) -> None:
        self.identity = identity


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
    """Parse scRNA Seq cell barcoding to assign cells to treatment, control or undetermined groups"""
    cell_identity_counts: dict = dict(
        Counter([cell.identity for cell in scrna_seq_data])
    )
    return cell_identity_counts


def secondary_demultiplexing(scrna_seq_data: list) -> dict:
    """
    Parse scRNA Seq cell barcoding to assign cells to treatment, control or undetermined groups
    after correcting for a Hamming Distance of 1
    """
    for cell in scrna_seq_data:
        # Iterate through cells that are undetermined & try to parse if their CellTags are actually control or treatment based upon minor Hamming Distance(s)
        control_barcode: str = cell.control_count[0]
        treatment_barcode: str = cell.treatment_count[0]
        if cell.identity == "undetermined":
            if cell.cell_tag1_count[1] == 1:
                cell_tag1_barcode: str = cell.cell_tag1_count[0]
                hamming_distance_treatment: int = calculate_hamming_distance(
                    treatment_barcode, cell_tag1_barcode
                )
                hamming_distance_control: int = calculate_hamming_distance(
                    control_barcode, cell_tag1_barcode
                )
                if hamming_distance_treatment <= 1:
                    cell.update_cell_identity("treatment")
                elif hamming_distance_control <= 1:
                    cell.update_cell_identity("control")
                else:
                    # Don't change cell identity from undetermined if it cannot be determined
                    continue
            elif cell.cell_tag2_count[1] == 1:
                cell_tag2_barcode: str = cell.cell_tag2_count[0]
                hamming_distance_treatment: int = calculate_hamming_distance(
                    treatment_barcode, cell_tag2_barcode
                )
                hamming_distance_control: int = calculate_hamming_distance(
                    control_barcode, cell_tag2_barcode
                )
                if hamming_distance_treatment <= 1:
                    cell.update_cell_identity("treatment")
                elif hamming_distance_control <= 1:
                    cell.update_cell_identity("control")
                else:
                    continue
            else:
                # Accounts for cases where no CellTags are present for a given Cell (0,0,0,0)
                cell.update_cell_identity("undetermined")
        else:
            continue
    cell_identity_counts = initial_demultiplexing(scrna_seq_data)
    return cell_identity_counts


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
secondary_cell_identity_counts: dict = secondary_demultiplexing(scrna_seq_data)
print("....Initial Demultiplexing....")
for cell_identity, count in cell_identity_counts.items():
    print(f"{cell_identity}: {count}")
print("....Secondary Demultiplexing....")
for cell_identity, count in secondary_cell_identity_counts.items():
    print(f"{cell_identity}: {count}")
