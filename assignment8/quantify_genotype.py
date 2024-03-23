from sys import argv, exit
from pathlib import Path
import itertools


SNV_FILE: Path = Path(argv[1])

"""
Python script to to quantify the number of homozygous and heterozygous SNVs and indels
for a given individual
python3 quantify_genotype.py <SNV VCF>
"""

if len(argv) != 2:
    print(__doc__)
    exit(1)


def SNV_VCF_IO(SNV_VCF: Path, individual_column: str):
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
        # Initialize different genotype conditions
        snv_genotypes_collection: dict = {}
        snv_genotypes_collection["homozygous_ref"] = 0
        snv_genotypes_collection["heterozygous"] = 0
        snv_genotypes_collection["homozygous_alt"] = 0
        snv_genotypes_collection["missing"] = 0
        for line in itertools.islice(vcf, header_start_index + 1, None):
            data: list = line.strip().split("\t")
            individual_data: list = data[individual_column_index].strip().split(":")
            genotype_data: str = individual_data[0]
            if genotype_data == "0/1":
                snv_genotypes_collection["heterozygous"] += 1
            elif genotype_data == "1/1":
                snv_genotypes_collection["homozygous_alt"] += 1
            elif genotype_data == "0/0":
                snv_genotypes_collection["homozygous_ref"] += 1
            else:
                snv_genotypes_collection["missing"] += 1
        return snv_genotypes_collection


def output_table_from_dict(data_set: dict) -> None:
    """Output a table from user provided dictionaries"""
    print("type\tcount")
    for key, value in data_set.items():
        print(f"{key}\t{value}")


SNV_GENOTYPES: dict = SNV_VCF_IO(SNV_FILE, "NA12878")
output_table_from_dict(SNV_GENOTYPES)
