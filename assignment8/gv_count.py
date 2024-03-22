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


def VCF_IO(VCF) -> pd.DataFrame:
    """Parse VCF & return dataframe of SV data"""
    # Parse raw file contents into list
    raw_file_contents: list = [
        line.strip().split("\t") for line in open(VCF).readlines()
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
    RAW_VCF_DF = pd.DataFrame(filtered_file_contents, columns=header_contents)
    return RAW_VCF_DF


SV_RAW_DF: pd.DataFrame = VCF_IO(SV_FILE)
SNV_RAW_DF: pd.DataFrame = VCF_IO(SNV_FILE)
