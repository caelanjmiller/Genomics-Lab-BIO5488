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
    FASTQ_FILE_EXTENSIONS: tuple = ('fastq', 'fq')
    FILENAME: str = os.path.basename(FILE)
    # Check to see if provided File argument is a valid FASTQ file
    if FILENAME.endswith(FASTQ_FILE_EXTENSIONS):
        # Create list of lines from FASTQ file & strip whitespace
        RAW_READS: list = [line.strip() for line in open(FILE).readlines()]
        # Split each FASTQ read into an element of a list --> every 4th line
        SPLIT_READS: list = [RAW_READS[line:line + 4] for line in range(0, len(RAW_READS), 4)]
        fastq: dict = {}
        # Iterate over each read & link identifier to both sequence & its quality score (as a dict)
        for READ in SPLIT_READS:
            fastq[READ[0]] = {READ[1].upper():READ[3]}
        return fastq
    else:
        raise Exception('Please provide a valid FAST file')
    

def count_nucleotides(fastq: dict) -> dict:
    """ Parse fastq as dict and return dict containing counts of canonical nucleotides """    
    # Initialize Counter subclass for tallying up nucleotides
    nucleotide_counter = Counter()
    # Iterate through FASTQ reads & update nucleotide count
    for read in fastq.values():
        for sequence in read.keys():
            nucleotide_counter.update(Counter(sequence))
    nucleotide_count = dict(nucleotide_counter)
    # Filter out non-canonical nucleotides from dinucleotide counts by iterating over list of dinucleotide keys ('AT', 'GC', etc) and removing them from dinucleotide count dictionary
    non_canonical_nucleotides: list = ['Y', 'N', 'W', 'R', 'B', 'H', 'V', 'S', 'K', 'D', 'M', 'U']
    for base in list(nucleotide_count.keys()):
        if any(nucleotide in base for nucleotide in non_canonical_nucleotides):
            del nucleotide_count[base]
    return dict(sorted(nucleotide_count.items()))


fastq = FASTQ_IO(FILE)
nucleotide_count = count_nucleotides(fastq)

if PRINTOUT == 'count':
    for nucleotide, count in nucleotide_count.items():
        print(f'{nucleotide}:{count}')