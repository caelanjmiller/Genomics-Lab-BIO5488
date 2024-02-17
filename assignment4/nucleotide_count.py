from sys import argv
import os
from collections import Counter
from itertools import islice
from statistics import mean 

FILE = argv[1]

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
        SEQUENCES: list = [line.strip().upper() for line in open(FILE).readlines()]
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


def count_dinucleotides(fasta: dict) -> dict:
    """Counts dinucleotides & returns a dict of counts"""
    fasta_dinucleotide_count: dict = {}
    for header, sequence in fasta.items():
        dinucleotides: list = []
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
        fasta_dinucleotide_count[header] = dinucleotide_count
    return fasta_dinucleotide_count


def calculate_dinucleotide_frequency(dinucleotide_count: dict) -> dict:
    """Calculate dinucleotide frequencies"""
    fasta_dinucleotide_frequencies: dict = {}
    for header, dinucleotide_count in dinucleotide_count.items():
        dinucleotide_frequency: dict = {}
        total_valid_dinucleotide_count: int = sum(dinucleotide_count.values())
        for dinucleotide, count in dinucleotide_count.items():
            dinucleotide_frequency[dinucleotide] = "{:.3f}".format(
                count / total_valid_dinucleotide_count
            )
        fasta_dinucleotide_frequencies[header] = dinucleotide_frequency
    return fasta_dinucleotide_frequencies


def printout_CpG_frequencies(dinucleotide_frequency: dict) -> float:
    """Print out CpG frequencies for multi FASTA"""
    cg_dinucleotide_frequencies: list = []
    for dinucleotide_frequencies_dict in dinucleotide_frequency.values():
        for dinucleotide, frequency in dinucleotide_frequencies_dict.items():
            if dinucleotide == "CG":
                cg_dinucleotide_frequencies.append(float(frequency))
    cg_dinucleotide_frequency: float = round(mean(cg_dinucleotide_frequencies), 3)
    return cg_dinucleotide_frequency

fasta: dict = FASTA_IO(FILE)
dinucleotide_count: dict = count_dinucleotides(fasta)
dinucleotide_frequency: dict = calculate_dinucleotide_frequency(dinucleotide_count)
cg_dinucleotide_frequency: float = printout_CpG_frequencies(dinucleotide_frequency)
print(cg_dinucleotide_frequency)

# Caelan Miller - 2024
