from sys import argv, exit
from pathlib import Path
import itertools
from collections import Counter
import matplotlib.pyplot as plt


SNV_FILE: Path = Path(argv[1])
SV_FILE: Path = Path(argv[2])

"""
Python script to count structural variants from provided VCFs
python3 gv_count.py <SNV VCF FILE> <SV VCF FILE>
"""


if len(argv) != 3:
    print(__doc__)
    exit(1)


def SV_VCF_IO(SV_VCF: Path, individual_column: str, data_needed: str) -> dict:
    """Parse SV VCF & return number of structural variants"""
    header_start_index: int = None
    header_contents: list = None
    individual_column_index: int = None

    with open(SV_VCF, "r") as vcf:
        # Find the index at which the header begins so we can skip forward in the file & collect header contents
        for index, line in enumerate(vcf):
            if "#CHROM" in line:
                header_start_index = index
                header_contents = line.strip().split("\t")
                individual_column_index = header_contents.index(individual_column)
                break
            else:
                continue

    with open(SV_VCF, "r") as vcf:
        if data_needed.upper() == "VARIANT":
            # Empty list for SV types to append to
            sv_types_collection: list = []
            for line in itertools.islice(vcf, header_start_index + 1, None):
                data: list = line.strip().split("\t")
                individual_data: list = data[individual_column_index].strip().split(":")
                if "0/1" in individual_data or "1/1" in individual_data:
                    sv_data_index: int = [
                        index for index, string in enumerate(data) if "SVTYPE" in string
                    ][0]
                    sv_data: list = data[sv_data_index].split(";")
                    sv_type_index: int = [
                        index
                        for index, string in enumerate(sv_data)
                        if "SVTYPE" in string
                    ][0]
                    sv_type: str = sv_data[sv_type_index].split("=")[1]
                    sv_types_collection.append(sv_type)
                else:
                    continue
            sv_types_counter: dict = dict(Counter(sv_types_collection))
            # 2 BND Calls represent 1 variant
            sv_types_counter["BND"] = int(sv_types_counter["BND"] / 2)
            return sv_types_counter
        elif data_needed.upper() == "LENGTH":
            # Empty list for SV lengths to append to
            sv_lengths_collection: dict = {}
            sv_lengths_collection["DEL"] = []
            sv_lengths_collection["MEI"] = []
            for line in itertools.islice(vcf, header_start_index + 1, None):
                data: list = line.strip().split("\t")
                individual_data: list = data[individual_column_index].strip().split(":")
                if "0/1" in individual_data or "1/1" in individual_data:
                    sv_data_index: int = [
                        index for index, string in enumerate(data) if "SVTYPE" in string
                    ][0]
                    sv_data: list = data[sv_data_index].split(";")
                    sv_type_index: int = [
                        index
                        for index, string in enumerate(sv_data)
                        if "SVTYPE" in string
                    ][0]
                    sv_type: str = sv_data[sv_type_index].split("=")[1]
                    if sv_type == "DEL":
                        sv_start_index: int = [
                            index
                            for index, string in enumerate(sv_data)
                            if "POS" in string
                        ][0]
                        sv_end_index: int = [
                            index
                            for index, string in enumerate(sv_data)
                            if "END" in string
                        ][0]
                        sv_start: int = int(sv_data[sv_start_index].split("POS=")[1])
                        sv_end: int = int(sv_data[sv_end_index].split("END=")[1])
                        sv_length: int = abs(sv_end - sv_start)
                        sv_lengths_collection["DEL"].append(sv_length)
                    elif sv_type == "MEI":
                        sv_start_index: int = [
                            index
                            for index, string in enumerate(sv_data)
                            if "POS" in string
                        ][0]
                        sv_end_index: int = [
                            index
                            for index, string in enumerate(sv_data)
                            if "END" in string
                        ][0]
                        sv_start: int = int(sv_data[sv_start_index].split("POS=")[1])
                        sv_end: int = int(sv_data[sv_end_index].split("END=")[1])
                        sv_length: int = abs(sv_end - sv_start)
                        sv_lengths_collection["MEI"].append(sv_length)
                else:
                    continue
            return sv_lengths_collection


def SNV_VCF_IO(SNV_VCF: Path, individual_column: str, data_needed: str):
    """Parse SNV VCF & return number of small genome variations (i.e. indels & SNPs)"""
    header_start_index: int = None
    header_contents: list = None
    individual_column_index: int = None

    with open(SNV_VCF, "r") as vcf:
        # Find the index at which the header begins so we can skip forward in the file & collect header contents
        for index, line in enumerate(vcf):
            if "#CHROM" in line:
                header_start_index = index
                header_contents = line.strip().split("\t")
                individual_column_index = header_contents.index(individual_column)
                break
            else:
                continue

    with open(SNV_VCF, "r") as vcf:
        if data_needed.upper() == "VARIANT":
            # Empty list for SNV types to append to
            snv_types_collection: list = []
            for line in itertools.islice(vcf, header_start_index + 1, None):
                data: list = line.strip().split("\t")
                individual_data: list = data[individual_column_index].strip().split(":")
                if "0/1" in individual_data or "1/1" in individual_data:
                    snv_ref_column_index: int = [
                        index
                        for index, string in enumerate(header_contents)
                        if "REF" in string
                    ][0]
                    snv_alt_column_index: int = [
                        index
                        for index, string in enumerate(header_contents)
                        if "ALT" in string
                    ][0]
                    alt_data: str = data[snv_alt_column_index]
                    ref_data: str = data[snv_ref_column_index]
                    # Calculate length of SNV or indel
                    snv_length: int = len(alt_data) - len(ref_data)
                    # Taking absolute value of SNV Length - signage determines biological relevance (neg indicates deletion, pos indicates insertion)
                    if abs(snv_length) >= 1:
                        snv_types_collection.append("INDEL")
                    else:
                        snv_types_collection.append("SNV")
            snv_types_counter: dict = dict(Counter(snv_types_collection))
            return snv_types_counter
        elif data_needed.upper() == "LENGTH":
            # Empty dict for SNV lengths to append to
            snv_lengths_collection: dict = {}
            snv_lengths_collection["INDEL"] = []
            for line in itertools.islice(vcf, header_start_index + 1, None):
                data: list = line.strip().split("\t")
                individual_data: list = data[individual_column_index].strip().split(":")
                if "0/1" in individual_data or "1/1" in individual_data:
                    snv_ref_column_index: int = [
                        index
                        for index, string in enumerate(header_contents)
                        if "REF" in string
                    ][0]
                    snv_alt_column_index: int = [
                        index
                        for index, string in enumerate(header_contents)
                        if "ALT" in string
                    ][0]
                    alt_data: str = data[snv_alt_column_index]
                    ref_data: str = data[snv_ref_column_index]
                    # Calculate absolute value of SNV or INDEL length & append to appropriate list
                    snv_length: int = abs(len(alt_data) - len(ref_data))
                    if abs(snv_length) >= 1:
                        snv_lengths_collection["INDEL"].append(snv_length)
                else:
                    continue
            return snv_lengths_collection


def output_table_from_dict(data_set1: dict, data_set2: dict) -> None:
    """Output a table from user provided dictionaries"""
    collated_data: dict = {**data_set1, **data_set2}
    print("type\tcount")
    for key, value in collated_data.items():
        print(f"{key}\t{value}")


def sv_proportion_of_total(sv_counter: dict, snv_counter: dict) -> float:
    """Calculate proportion of genomic variation that are SVs"""
    collated_data: dict = {**sv_counter, **snv_counter}
    total_count: int = sum(list(collated_data.values()))
    sv_proportion: float = round((sum(list(sv_counter.values())) / total_count), 4)
    return sv_proportion


def generate_histogram_sv_genomic_variation(
    sv_data: dict, title: tuple, filename: str
) -> None:
    """Generate histogram of lengths for genomic variations"""
    current_directory: Path = Path.cwd()
    plt.figure()
    plt.hist(list(sv_data.values()), bins=20)
    plt.autoscale()
    plt.xlabel("Length (bp)")
    plt.ylabel("Count")
    plt.title(f"Length of {title[0]}s in {title[1]}")
    plt.savefig(f"{current_directory}/histogram_{filename}.png")
    plt.close()


SV_COUNTER: dict = SV_VCF_IO(SV_FILE, "NA12878", "variant")
SNV_COUNTER: dict = SNV_VCF_IO(SNV_FILE, "NA12878", "variant")
output_table_from_dict(SV_COUNTER, SNV_COUNTER)
sv_proportion: float = sv_proportion_of_total(SV_COUNTER, SNV_COUNTER)
print(sv_proportion)

SV_LENGTH_COUNTER: dict = SV_VCF_IO(SV_FILE, "NA12878", "length")
SNV_LENGTH_COUNTER: dict = SNV_VCF_IO(SNV_FILE, "NA12878", "length")
SV_LENGTHS: list = [{sv: lengths} for sv, lengths in SV_LENGTH_COUNTER.items()]
generate_histogram_sv_genomic_variation(
    SV_LENGTHS[0], tuple(("DEL", "NA12878")), "deletions"
)
generate_histogram_sv_genomic_variation(
    SV_LENGTHS[1], tuple(("MEI", "NA12878")), "meis"
)
generate_histogram_sv_genomic_variation(
    SNV_LENGTH_COUNTER, tuple(("INDEL", "NA12878")), "indels"
)
