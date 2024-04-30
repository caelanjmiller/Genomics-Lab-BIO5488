#!/usr/bin/env python3

"""
Python script to calculate the fraction of wobble positions that
are conserved across all species in a user provided CLUSTAL W alignment file
Usage: python3 neutral_rate.py <clustalw alignment file>
"""

import os
from itertools import islice
from sys import argv
from pathlib import Path

CLUSTALW_ALN: Path = Path(argv[1])
FASTA_SEQ: Path = Path(argv[2])


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
            fasta_header = FASTA_ENTRY[0].replace(">", "")
            fasta_sequence = "".join(FASTA_ENTRY[1:])
            fasta[fasta_header] = fasta_sequence
        return fasta
    else:
        raise Exception("Please provide a valid FASTA file")


def extract_CDS_sequence(fasta_sequences: dict) -> dict:
    """Extract CDS Sequences from genomic sequences"""
    cds_sequences: dict = {}
    for header, genomic_sequence in fasta_sequences.items():
        cds_start = None
        for index, nucleotide in enumerate(list(genomic_sequence)):
            if nucleotide.isupper():
                cds_start = index
                # Break at 1st occurrence of capital letter - A in ATG
                break
        cds_sequences[header] = genomic_sequence[cds_start:]
    return cds_sequences


fasta_sequences: dict = FASTA_IO(FASTA_SEQ)
cds_sequences: dict = extract_CDS_sequence(fasta_sequences)
for name, sequence in cds_sequences.items():
    print(f"{name}: {len(sequence)}")
