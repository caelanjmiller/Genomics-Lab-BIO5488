from sys import argv
import os
import math
import re
from collections import Counter
import matplotlib.pyplot as plt
from pathlib import Path

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
        fasta[SEQUENCE[0].replace(">", "")] = "".join(SEQUENCE[1:]).upper()
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
        nucleotide_frequency[nucleotide] = "{:.3f}".format(
            count / total_valid_dinucleotide_count
        )
    return dict(sorted(nucleotide_frequency.items()))


def calculate_dinucleotide_frequency(dinucleotide_count: dict) -> dict:
    """Calculate dinucleotide frequencies"""
    dinucleotide_frequency: dict = {}
    total_valid_dinucleotide_count: int = sum(dinucleotide_count.values())
    for dinucleotide, count in dinucleotide_count.items():
        dinucleotide_frequency[dinucleotide] = "{:.3f}".format(
            count / total_valid_dinucleotide_count
        )
    return dict(sorted(dinucleotide_frequency.items()))


def find_gaps(fasta: dict) -> list:
    """Find gaps (defined by N) in a provided FASTA sequence and return a list of each gapped sequence (string of Ns)"""
    n_gaps: list = []
    # Iterate through FASTA dict for sequence string
    for header, sequence in fasta.items():
        # Use regular expressions to find occurrences of N & create a tuple of start index & stop index
        for gap in re.finditer(r"N+", sequence):
            n_gap_indices: tuple = tuple((gap.start(), gap.end()))
            n_gaps.append(n_gap_indices)
    return n_gaps


def find_gap_lengths(gaps: list) -> list:
    """Calculate sequence length of each gap"""
    gap_lengths: list = []
    for gap in gaps:
        start, stop = gap
        gap_length: int = int(stop - start)
        gap_lengths.append(gap_length)
    return gap_lengths


def create_gap_length_histogram(gap_lengths: list, FILE):
    """Create basic histogram of gap length for a given FASTA"""
    FILENAME: str = os.path.basename(FILE).split(".fasta")[0]
    current_directory: Path = Path.cwd()
    log10_gap_lengths: list = [math.log10(gap) for gap in gap_lengths]
    plt.hist(log10_gap_lengths, bins=20)
    plt.ylabel("Frequency")
    plt.xlabel("Gap Length (bp) [log10]")
    plt.title(f"Gap Lengths in {FILENAME} Assembly")
    plt.savefig(f"{current_directory}/{FILENAME}_gap_distribution.png")


fasta = FASTA_IO(FILE)
nucleotide_count = count_nucleotides(fasta)
nucleotide_frequency = calculate_nucleotide_frequency(nucleotide_count)
dinucleotide_count = count_dinucleotides(fasta)
dinucleotide_frequency = calculate_dinucleotide_frequency(dinucleotide_count)
genome_length = calculate_sequence_length(fasta)
gaps = find_gaps(fasta)
gap_lengths = find_gap_lengths(gaps)

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
elif PRINTOUT == "length":
    for length in genome_length.values():
        print(f"Genome is {length} bp")
elif PRINTOUT == "graph":
    create_gap_length_histogram(gap_lengths, FILE)
    print(len(gap_lengths))
    gap_counter: dict = dict(Counter(gap_lengths))
    print("Length(bp): Count")
    for length, count in gap_counter.items():
        print(f"{length}: {count}")
else:
    raise Exception("Provide valid printout option")
