from sys import argv
import os
from collections import Counter

FILE = argv[1]

def FASTQ_IO(FILE) -> dict:
    """ Parse in valid FASTQ file & return d"""
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
            fastq[READ[0]] = {READ[1]:READ[3]}
        return fastq
    else:
        raise Exception('Please provide a valid FAST file')
    

fastq = FASTQ_IO(FILE)