from sys import argv
import os
from collections import Counter
import seaborn as sns

FILE = argv[1]
PRINTOUT = argv[2]

def FASTA_IO(FILE) -> dict:
    """Parse in valid FASTA file & return dictionary containing header as key and str of sequence as value"""
    # Valid FASTA Extensions
    FASTA_FILE_EXTENSIONS: tuple = (
        "fasta",
        "fas",
        "fa",
        "fna",
        "ffn",
        "faa",
        "mpfa",
        "frn",
    )
    FILENAME: str = os.path.basename(FILE)
    # Check to see if provided File argument is a valid FASTA file
    if FILENAME.endswith(FASTA_FILE_EXTENSIONS):
        # Create list of lines from FASTA file & strip whitespace
        SEQUENCE: list = [line.strip() for line in open(FILE).readlines()]
        fasta: dict = {}
        fasta[SEQUENCE[0].replace(">", "")] = "".join(SEQUENCE[1:])
        return fasta
    else:
        raise Exception("Please provide a valid FASTA file")