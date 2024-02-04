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


def calculate_sequence_length(fasta: dict) -> dict:
    """Calculate length of a given FASTA sequence"""
    sequence_length: dict = {}
    for header, sequence in fasta.items():
        sequence_length[header] = len(sequence)
    return sequence_length


def count_nucleotides(fasta: dict) -> dict:
    """Parse FASTA as dict and return dict containing counts of canonical nucleotides"""
    sequence_str: str = ""
    for sequence in fasta.values():
        sequence_str += sequence
    # Utilize Counter subclass for counting nucleotides
    nucleotide_counter = Counter(sequence_str)
    nucleotide_count: dict = dict(nucleotide_counter)
    # Filter out non-canonical nucleotides from nucleotide counts by iterating over list of nucleotide keys
    # and removing them from nucleotide count dictionary
    non_canonical_nucleotides: list = [
        "Y",
        "N",
        "W",
        "R",
        "B",
        "H",
        "V",
        "S",
        "K",
        "D",
        "M",
        "U",
    ]
    for base in list(nucleotide_count.keys()):
        if any(nucleotide in base for nucleotide in non_canonical_nucleotides):
            del nucleotide_count[base]
    return dict(sorted(nucleotide_count.items()))


def count_dinucleotides(fasta: dict) -> dict:
    """Counts dinucleotides & returns a dict of counts"""
    dinucleotides: list = []
    fasta_header: str = list(fasta.keys())[0]
    # Convert sequence str from fasta dict into a list of individual characters
    sequence: list = [*fasta[fasta_header]]
    # Iterate through sequence & take slice of sequence list (2 nucleotides) and join together as a string and append to list
    for nucleotide_index in range(len(sequence) - 1):
        dinucleotide_slice: str = "".join(
            sequence[nucleotide_index : (nucleotide_index + 2)]
        )
        dinucleotides.append(dinucleotide_slice)
    # Utilize Counter subclass to tally up the dinucleotides
    dinucleotide_count: dict = dict(Counter(dinucleotides))
    # Filter out non-canonical nucleotides from dinucleotide counts by iterating over list of dinucleotide keys ('AT', 'GC', etc) and removing them from dinucleotide count dictionary
    non_canonical_nucleotides: list = [
        "Y",
        "N",
        "W",
        "R",
        "B",
        "H",
        "V",
        "S",
        "K",
        "D",
        "M",
        "U",
    ]
    for dinucleotide in list(dinucleotide_count.keys()):
        if any(nucleotide in dinucleotide for nucleotide in non_canonical_nucleotides):
            del dinucleotide_count[dinucleotide]
    return dict(sorted(dinucleotide_count.items()))


def calculate_nucleotide_frequency(nucleotide_count: dict) -> dict:
    """Calculate nucleotide frequencies"""
    nucleotide_frequency: dict = {}
    total_valid_dinucleotide_count: int = sum(nucleotide_count.values())
    for nucleotide, count in nucleotide_count.items():
        nucleotide_frequency[nucleotide]: str = "{:.3f}".format(
            count / total_valid_dinucleotide_count
        )
    return dict(sorted(nucleotide_frequency.items()))


def calculate_dinucleotide_frequency(dinucleotide_count: dict) -> dict:
    """Calculate dinucleotide frequencies"""
    dinucleotide_frequency: dict = {}
    total_valid_dinucleotide_count: int = sum(dinucleotide_count.values())
    for dinucleotide, count in dinucleotide_count.items():
        dinucleotide_frequency[dinucleotide]: str = "{:.3f}".format(
            count / total_valid_dinucleotide_count
        )
    return dict(sorted(dinucleotide_frequency.items()))
