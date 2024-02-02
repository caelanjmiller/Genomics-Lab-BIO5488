from sys import argv
import os
from collections import Counter

FILE = argv[1]
PRINTOUT = argv[2]


def FASTQ_IO(FILE) -> dict:
    """
    Parse in valid FASTQ file & return dictionary of identifier with sequence & its quality scores
    """
    # Valid FASTQ Extensions
    FASTQ_FILE_EXTENSIONS: tuple = ("fastq", "fq")
    FILENAME: str = os.path.basename(FILE)
    # Check to see if provided File argument is a valid FASTQ file
    if FILENAME.endswith(FASTQ_FILE_EXTENSIONS):
        # Create list of lines from FASTQ file & strip whitespace
        RAW_READS: list = [line.strip() for line in open(FILE).readlines()]
        # Split each FASTQ read into an element of a list --> every 4th line
        SPLIT_READS: list = [
            RAW_READS[line : line + 4] for line in range(0, len(RAW_READS), 4)
        ]
        fastq: dict = {}
        # Iterate over each read & link identifier to both sequence & its quality score (as a dict)
        for READ in SPLIT_READS:
            fastq[READ[0]] = {READ[1].upper(): READ[3]}
        return fastq
    else:
        raise Exception("Please provide a valid FAST file")


def count_nucleotides(fastq: dict) -> dict:
    """Parse fastq as dict and return dict containing counts of canonical nucleotides"""
    # Initialize Counter subclass for tallying up nucleotides
    nucleotide_counter = Counter()
    # Iterate through FASTQ reads & update nucleotide count
    for read in fastq.values():
        for sequence in read.keys():
            nucleotide_counter.update(Counter(sequence))
    nucleotide_count = dict(nucleotide_counter)
    # Filter out non-canonical nucleotides from nucleotide counts by iterating over list of nucleotide keys and removing them from nucleotide count dictionary
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


def count_dinucleotides(fastq: dict) -> dict:
    """Counts dinucleotides & returns a dict of counts"""
    dinucleotides: list = []
    for read in fastq.values():
        for sequence in read.keys():
            for nucleotide_index in range(len(sequence) - 1):
                dinucleotide_slice: str = "".join(
                    sequence[nucleotide_index : (nucleotide_index + 2)]
                )
                dinucleotides.append(dinucleotide_slice)
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
    total_valid_basecount: int = sum(nucleotide_count.values())
    nucleotide_frequencies: dict = {}
    for nucleotide, count in nucleotide_count.items():
        nucleotide_frequencies[nucleotide] = round((count / total_valid_basecount), 3)
    return nucleotide_frequencies


def calculate_dinucleotide_frequency(dinucleotide_count: dict) -> dict:
    """Calculate dinucleotide frequencies"""
    dinucleotide_frequency: dict = {}
    total_valid_dinucleotide_count: int = sum(dinucleotide_count.values())
    for dinucleotide, count in dinucleotide_count.items():
        dinucleotide_frequency[dinucleotide] = round(
            (count / total_valid_dinucleotide_count), 3
        )
    return dinucleotide_frequency


fastq = FASTQ_IO(FILE)
nucleotide_count = count_nucleotides(fastq)
dinucleotide_count = count_dinucleotides(fastq)
nucleotide_frequency = calculate_nucleotide_frequency(nucleotide_count)
dinucleotide_frequency = calculate_dinucleotide_frequency(dinucleotide_count)

if PRINTOUT == "count":
    for nucleotide, count in nucleotide_count.items():
        print(f"{nucleotide}:{count}")
elif PRINTOUT == "dicount":
    for dinucleotide, count in dinucleotide_count.items():
        print(f"{dinucleotide}:{count}")
elif PRINTOUT == "frequency":
    for nucleotide, frequency in nucleotide_frequency.items():
        print(f"{nucleotide}:{frequency}")
elif PRINTOUT == "difrequency":
    for dinucleotide, frequency in dinucleotide_frequency.items():
        print(f"{dinucleotide}:{frequency}")
else:
    raise Exception("Provide a valid printout option")

# Caelan Miller - 2024
