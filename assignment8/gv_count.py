from sys import argv, exit
from pathlib import Path
import pandas as pd

SNV_FILE: Path = Path(argv[1])
SV_FILE: Path = Path(argv[2])

"""
python3 gv_count.py <SNV VCF FILE> <SV VCF FILE>
"""


if len(argv) != 3:
    print(__doc__)
    exit(1)


def SV_IO(SV_VCF) -> pd.DataFrame:
    """Parse SV VCF & return dataframe of SV data"""
    # Parse raw file contents into list
    raw_file_contents: list = [
        line.strip().split("\t") for line in open(SV_VCF).readlines()
    ]
    # Get index of header
    header_start_index: int = [
        index for index, line in enumerate(raw_file_contents) if "#CHROM" in line
    ][0]
    # Create list of column names
    header_contents: list = raw_file_contents[header_start_index]
    # Filter (skip) rows of data in VCF & create list of data
    filtered_file_contents: list = raw_file_contents[header_start_index + 1 :]
    # Initialize a dataframe of the data and return it
    RAW_SV_DF = pd.DataFrame(filtered_file_contents, columns=header_contents)
    return RAW_SV_DF


def SNV_IO(SNV_VCF) -> list:
    """Parse SNV VCF & return dataframe of SNV data"""
    pass


SV_RAW_DF: pd.DataFrame = SV_IO(SV_FILE)
