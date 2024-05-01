#!/usr/bin/env python3

"""
Python script to calculate the fraction of wobble positions that
are conserved across all species in a user provided CLUSTAL W alignment file
Usage: python3 neutral_rate.py <clustalw alignment file> <FASTA file>
"""

import scipy
import itertools
from sys import argv
from pathlib import Path
from Bio import SeqIO, AlignIO
from textwrap import wrap

import scipy.stats

CLUSTALW_ALN: Path = Path(argv[1])
FASTA_SEQ: Path = Path(argv[2])


def FASTA_IO(FILE) -> dict:
    """Parse in valid FASTA file & return dictionary containing header as key and str of sequence as value"""
    fasta_sequences: dict = {}
    for sequence_record in SeqIO.parse(FILE, "fasta"):
        fasta_sequences[sequence_record.id] = str(sequence_record.seq)
    return fasta_sequences


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


def alignment_IO(alignment_file: Path):
    """Parse CLUSTAL alignment file & return an BioPython alignment object"""
    alignment = AlignIO.read(alignment_file, "clustal")
    return alignment


def nucleotide_split_codons(cds_sequences: dict) -> dict:
    """Split CDS Sequences into putative codons (sets of 3)"""
    codon_sequences: dict = {}
    for header, sequence in cds_sequences.items():
        codons = wrap(sequence, 3)
        codon_sequences[header] = codons
    return codon_sequences


def calculate_wobble_conservation(codon_sequences: dict) -> float:
    """Calculate fraction of conserved wobble positions across CDS sequences"""
    third_nucleotide_str_dict: dict = {}
    for header, codon_sequence in codon_sequences.items():
        third_string: list = []
        for codon in codon_sequence:
            third_string.append(codon[2])
        third_nucleotide_str_dict[header] = third_string
    total_wobble_positions: int = int(
        (len(next(iter(third_nucleotide_str_dict.values()))))
    )
    sequence_one: list = third_nucleotide_str_dict["Scer"]
    sequence_two: list = third_nucleotide_str_dict["Skud"]
    sequence_three: list = third_nucleotide_str_dict["Smik"]
    sequence_four: list = third_nucleotide_str_dict["Sbay"]
    conserved_positions: int = 0
    for i, j, x, z in itertools.zip_longest(
        sequence_one, sequence_two, sequence_three, sequence_four
    ):
        if i == j and i == x and i == z:
            conserved_positions += 1
        else:
            continue
    wobble_conservation: float = round(
        float(conserved_positions / total_wobble_positions), 4
    )
    return wobble_conservation


CLUSTALW_alignment = alignment_IO(CLUSTALW_ALN)
fasta_sequences: dict = FASTA_IO(FASTA_SEQ)
cds_sequences: dict = extract_CDS_sequence(fasta_sequences)
codon_sequences: dict = nucleotide_split_codons(cds_sequences)
wobble_conservation: float = calculate_wobble_conservation(codon_sequences)
print(f"Conserved wobble positions: {wobble_conservation}")
print(
    f"Number of bp for a 10bp sequence to be more conserved than expected: {int(scipy.stats.binom.ppf(.95, 10, wobble_conservation))}"
)
