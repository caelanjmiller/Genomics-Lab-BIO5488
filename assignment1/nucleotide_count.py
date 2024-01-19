from sys import argv
import os

FILE = argv[1]

def FASTA_IO(FILE) -> dict:
    """ Parse in valid FASTA file & return dictionary containing header as key and str of sequence as value """
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


def nucleotide_count(fasta: dict) -> dict:
    """ Counts nucleotides & returns a dict of counts """
    base_count: dict = {'A':0, 'T':0, 'G': 0, 'C': 0}
    fasta_header: str = list(fasta.keys())[0]
    # Convert sequence str from fasta dict into a list of individual characters
    sequence: list = [*fasta[fasta_header]]
    for nucleotide in sequence:
        match nucleotide:
            case 'A' | 'a':
                base_count['A'] += 1
            case 'T' | 't':
                base_count['T'] += 1
            case 'G' | 'g':
                base_count['G'] += 1
            case 'C' | 'c':
                base_count['C'] += 1
            case _:
                continue
    return base_count

fasta = FASTA_IO(FILE)
print(nucleotide_count(fasta))