from sys import argv
from pathlib import Path

FREQUENCY_FILE1 = Path(argv[1])
FREQUENCY_FILE2 = Path(argv[2])


def write_file_contents_to_dict(file: Path) -> dict:
    """Write file contents to dict for downstream analysis"""
    FILE_CONTENTS: list = [line.strip() for line in open(file).readlines()]
    freq_data: dict = {}
    # Iterate through list of file contents, assigning the data to a dictionary
    for line in FILE_CONTENTS:
        data_split: list = line.split(":")
        key: str = data_split[0]
        value: float = float(data_split[1])
        freq_data[key] = value
    return freq_data


def enrichment_calculations(freq_one: dict, freq_two: dict) -> dict:
    """Calculate enrichment of the values from both dicts"""
    enrichment_score: dict = {}
    for key_one, frequency_one in freq_one.items():
        for key_two, frequency_two in freq_two.items():
            if key_one == key_two:
                enrichment_score[key_one] = round(
                    float(frequency_one / frequency_two), 3
                )
    return enrichment_score


frequency_one: dict = write_file_contents_to_dict(FREQUENCY_FILE1)
frequency_two: dict = write_file_contents_to_dict(FREQUENCY_FILE2)
enrichment_scores: dict = enrichment_calculations(frequency_one, frequency_two)

for key, frequency in enrichment_scores.items():
    print(f"{key}:{frequency}")

# Caelan Miller - 2024
