from sys import argv
import os

FILE = argv[1]
PRINTOUT = argv[2]

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

def count_nucleotides(fasta: dict) -> dict:
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

def calculate_nucleotide_frequency(base_count: dict) -> dict:
    """Calculate nucleotide frequencies from base counts"""
    total_valid_basecount: int = sum(base_count.values())
    nucleotide_frequencies: dict = {}
    # Calculate nucleotide frequency and format to 3 decimal places
    nucleotide_frequencies['A']: str = "{:.3f}".format(base_count['A'] / total_valid_basecount)
    nucleotide_frequencies['T']: str = "{:.3f}".format(base_count['T'] / total_valid_basecount)
    nucleotide_frequencies['G']: str = "{:.3f}".format(base_count['G'] / total_valid_basecount)
    nucleotide_frequencies['C']: str = "{:.3f}".format(base_count['C'] / total_valid_basecount)
    return nucleotide_frequencies

fasta = FASTA_IO(FILE)
nucleotide_count = count_nucleotides(fasta)
nucleotide_frequency = calculate_nucleotide_frequency(nucleotide_count)

if PRINTOUT == 'count':
    print(nucleotide_count)
elif PRINTOUT == 'frequency':
    print(nucleotide_frequency)
else:
    raise Exception('Provide valid printout option')