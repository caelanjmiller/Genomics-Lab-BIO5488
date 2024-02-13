from sys import argv
import os
from collections import Counter
from itertools import islice

FILE = argv[1]
FASTA_HEADER = argv[2]
PRINTOUT = argv[3]


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
        SEQUENCES: list = [line.strip() for line in open(FILE).readlines()]
        fasta: dict = {}
        # Capture indices at which a FASTA header is seen
        fasta_header_indices: list = [
            index for index, line in enumerate(SEQUENCES) if line.startswith(">")
        ]
        # Capture FASTA Header & Sequence as a list of lists
        sequence_information_indices: list = [
            tuple((fasta_header_indices[index], fasta_header_indices[index + 1]))
            for index in range(len(fasta_header_indices) - 1)
        ]
        captured_fasta_information: list = [
            list(islice(SEQUENCES, index, sequences))
            for index, sequences in sequence_information_indices
        ]
        for FASTA_ENTRY in captured_fasta_information:
            fasta_header = FASTA_ENTRY[0]
            fasta_sequence = "".join(FASTA_ENTRY[1:])
            fasta[fasta_header] = fasta_sequence
        return fasta
    else:
        raise Exception("Please provide a valid FASTA file")


def create_appended_fasta_sequence(fasta: dict, header: str) -> dict:
    """Create large appended FASTA sequence"""
    mega_fasta: dict = {}
    sequence_str: str = ""
    for sequence in fasta.values():
        sequence_str += sequence
    mega_fasta[header] = sequence_str
    return mega_fasta


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


def calculate_nucleotide_frequency(base_count: dict) -> dict:
    """Calculate nucleotide frequencies from base counts"""
    total_valid_basecount: int = sum(base_count.values())
    nucleotide_frequencies: dict = {}
    # Calculate nucleotide frequency and format to 3 decimal places
    nucleotide_frequencies["A"]: str = "{:.3f}".format(
        base_count["A"] / total_valid_basecount
    )
    nucleotide_frequencies["T"]: str = "{:.3f}".format(
        base_count["T"] / total_valid_basecount
    )
    nucleotide_frequencies["G"]: str = "{:.3f}".format(
        base_count["G"] / total_valid_basecount
    )
    nucleotide_frequencies["C"]: str = "{:.3f}".format(
        base_count["C"] / total_valid_basecount
    )
    return dict(sorted(nucleotide_frequencies.items()))


def calculate_dinucleotide_frequency(dinucleotide_count: dict) -> dict:
    """Calculate dinucleotide frequencies"""
    dinucleotide_frequency: dict = {}
    total_valid_dinucleotide_count: int = sum(dinucleotide_count.values())
    for dinucleotide, count in dinucleotide_count.items():
        dinucleotide_frequency[dinucleotide]: str = "{:.3f}".format(
            count / total_valid_dinucleotide_count
        )
    return dict(sorted(dinucleotide_frequency.items()))


fasta = FASTA_IO(FILE)
mega_fasta = create_appended_fasta_sequence(fasta, FASTA_HEADER)
nucleotide_count = count_nucleotides(mega_fasta)
nucleotide_frequency = calculate_nucleotide_frequency(nucleotide_count)
dinucleotide_count = count_dinucleotides(mega_fasta)
dinucleotide_frequency = calculate_dinucleotide_frequency(dinucleotide_count)

if PRINTOUT == "count":
    for nucleotide, count in nucleotide_count.items():
        print(f"{nucleotide}:{count}")
elif PRINTOUT == "frequency":
    for nucleotide, frequency in nucleotide_frequency.items():
        print(f"{nucleotide}:{frequency}")
elif PRINTOUT == "dicount":
    for dinucleotide, count in dinucleotide_count.items():
        print(f"{dinucleotide}:{count}")
elif PRINTOUT == "difrequency":
    for dinucleotide, frequency in dinucleotide_frequency.items():
        print(f"{dinucleotide}:{frequency}")
else:
    raise Exception("Provide valid printout option")

# Caelan Miller - 2024
