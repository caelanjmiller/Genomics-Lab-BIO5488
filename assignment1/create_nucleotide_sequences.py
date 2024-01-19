from sys import argv

LENGTH: int = argv[1]
A_FREQ: float = argv[2]
C_FREQ: float = argv[3]
G_FREQ: float = argv[4]
T_FREQ: float = argv[5]
FASTA_OUTPUT: str = argv[6]

NUCLEOTIDE_FREQ: dict = {'A': A_FREQ, 'T': T_FREQ, 'C': C_FREQ, 'G': G_FREQ}

def create_novel_dna_sequence(length: int, frequencies: dict) -> dict:
    """ 
    Create a novel DNA sequence that is of a provided length & 
    similar to provided nucleotide frequencies (within 2 decimal places)
    """
    generated_dna_sequence: str = ''
    novel_seq: dict = {}
    novel_seq['GENERATED_SEQUENCE'] = generated_dna_sequence
    return novel_seq

def sequence_to_FASTA (novel_seq: dict, FASTA_OUTPUT: str):
    """
    Convert novel generated DNA sequence to FASTA format
    """
    fasta_header: str = list(novel_seq.keys())[0]
    with open(f'{FASTA_OUTPUT}.fa', 'w') as FASTA:
        FASTA.write(f'>{fasta_header}\n')
        FASTA.write(novel_seq[fasta_header])


dna_sequence = create_novel_dna_sequence(LENGTH, NUCLEOTIDE_FREQ)
sequence_to_FASTA(dna_sequence, FASTA_OUTPUT)