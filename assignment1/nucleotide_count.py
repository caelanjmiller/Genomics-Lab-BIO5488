from sys import argv
import os

FILE = argv[1]

def FASTA_IO(FILE) -> dict:
    """ Parse in valid FASTA file & return list of sequence """
    # Valid FASTA Extensions
    FASTA_FILE_EXTENSIONS: tuple = ('fasta', 'fas', 'fa', 'fna', 'ffn','faa','mpfa','frn')
    FILENAME: str = os.path.basename(FILE)
    # Check to see if provided File argument is a valid FASTA file
    if FILENAME.endswith(FASTA_FILE_EXTENSIONS):
        SEQUENCE: list = [line.strip() for line in open(FILE).readlines()]
        fasta: dict = {}
        fasta[SEQUENCE[0].replace('>', '')] = "".join(SEQUENCE[1:])
        return fasta
    else:
        raise Exception('Please provide a valid FASTA file')


def nucleotide_count (SEQUENCE: list) -> dict:
    """ Calculates nucleotide frequency & returns a dict of frequencies """
    pass
