from sys import argv
import os
from collections import Counter

FILE = argv[1]
PRINTOUT = argv[2]

def FASTA_IO(FILE) -> dict:
    """ Parse in valid FASTA file & return dictionary containing header as key and str of sequence as value """
    # Valid FASTA Extensions
    FASTA_FILE_EXTENSIONS: tuple = ('fasta', 'fas', 'fa', 'fna', 'ffn','faa','mpfa','frn')
    FILENAME: str = os.path.basename(FILE)
    # Check to see if provided File argument is a valid FASTA file
    if FILENAME.endswith(FASTA_FILE_EXTENSIONS):
        # Create list of lines from FASTA file & strip whitespace
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

def count_dinucleotides(fasta: dict) -> dict:
    """ Counts dinucleotides & returns a dict of counts """
    dinucleotides: list = []
    fasta_header: str = list(fasta.keys())[0]
    # Convert sequence str from fasta dict into a list of individual characters
    sequence: list = [*fasta[fasta_header]]
    # Iterate through sequence & take slice of sequence list (2 nucleotides) and join together as a string and append to list
    for nucleotide_index in range(len(sequence) - 1):
        dinucleotide_slice: str = "".join(sequence[nucleotide_index:(nucleotide_index + 2)])
        dinucleotides.append(dinucleotide_slice)
    # Utilize Counter subclass to tally up the dinucleotides
    dinucleotide_count: dict = dict(Counter(dinucleotides))
    return dinucleotide_count

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

def calculate_dinucleotide_frequency(base_count: dict) -> dict:
    """ Calculate dinucleotide frequencies"""
    pass


fasta = FASTA_IO(FILE)
nucleotide_count = count_nucleotides(fasta)
nucleotide_frequency = calculate_nucleotide_frequency(nucleotide_count)
dinucleotide_count = count_dinucleotides(fasta)
dinucleotide_frequency = calculate_dinucleotide_frequency(dinucleotide_count)

if PRINTOUT == 'count':
    print(nucleotide_count)
elif PRINTOUT == 'frequency':
    print(nucleotide_frequency)
elif PRINTOUT == 'dicount':
    print(dinucleotide_count)
elif PRINTOUT == 'difrequency':
    print(dinucleotide_frequency)
else:
    raise Exception('Provide valid printout option')