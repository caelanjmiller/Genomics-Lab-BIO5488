from pathlib import Path
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


def parse_cell_matrix_csv(FILE: Path) -> list:
    """Parse cell matrix CSV & initialize Cell objects from provided information"""
    with open(FILE, "r") as CSV:
        scrna_seq_data: list = []
        file_contents: list = csv.reader(CSV, delimiter=",")
        csv_header: list = next(file_contents)
        for line in file_contents:
            cell: Cell = Cell(
                line[0],
                tuple(line[1]),
                tuple((csv_header[1], line[2])),
                tuple((csv_header[2], line[3])),
                tuple((csv_header[3], line[4])),
            )
            scrna_seq_data.append(cell)
    return scrna_seq_data


def initial_demultiplexing(scrna_seq_data: list) -> dict:
    """Parse scRNA Seq cell barcoding to assign cells to treatment, control or nondetermined groups"""
    pass


scrna_seq_data: list = parse_cell_matrix_csv(CALL_MATRIX)
