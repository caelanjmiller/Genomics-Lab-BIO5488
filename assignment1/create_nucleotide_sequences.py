from sys import argv
import random

# Parse user input & initialize variable to correct data type
LENGTH: int = int(argv[1])
A_FREQ: float = float(argv[2])
C_FREQ: float = float(argv[3])
G_FREQ: float = float(argv[4])
T_FREQ: float = float(argv[5])
FASTA_OUTPUT: str = argv[6]

NUCLEOTIDE_FREQ: dict = {"A": A_FREQ, "C": C_FREQ, "G": G_FREQ, "T": T_FREQ}


def create_novel_dna_sequence(length: int, frequencies: dict) -> dict:
    """
    Create a novel DNA sequence that is of a provided length &
    similar to provided nucleotide frequencies (within 2 decimal places)
    """
    total_nucleotide_frequency = sum(frequencies.values())
    IDEAL_NUCLEOTIDE_FREQ: float = float(1.0)
    if (
        total_nucleotide_frequency > IDEAL_NUCLEOTIDE_FREQ
        or total_nucleotide_frequency < IDEAL_NUCLEOTIDE_FREQ
    ):
        raise Exception("Provided nucleotide frequencies do not sum to 1")
    else:
        # Creates a DNA sequence string of a provided length that from the four choices are weighted (to be chosen at random) by the provided frequencies
        generated_dna_sequence: str = "".join(
            random.choices(
                population=["A", "C", "G", "T"], weights=frequencies.values(), k=length
            )
        )
        novel_seq: dict = {}
        novel_seq["GENERATED_SEQUENCE"] = generated_dna_sequence
        return novel_seq


def sequence_to_FASTA(novel_seq: dict, FASTA_OUTPUT: str):
    """
    Convert novel generated DNA sequence to FASTA format
    """
    fasta_header: str = list(novel_seq.keys())[0]
    with open(f"{FASTA_OUTPUT}.fa", "w") as FASTA:
        FASTA.write(f">{fasta_header}\n")
        FASTA.write(novel_seq[fasta_header])


dna_sequence = create_novel_dna_sequence(LENGTH, NUCLEOTIDE_FREQ)
sequence_to_FASTA(dna_sequence, FASTA_OUTPUT)
