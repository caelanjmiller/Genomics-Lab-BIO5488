#!/usr/bin/env python3

"""
Python script to calculate the fraction of wobble positions that
are conserved across all species in a user provided CLUSTAL W alignment file
Usage: python3 neutral_rate.py <clustalw alignment file>
"""

from sys import argv
from pathlib import Path

CLUSTALW_ALN: Path = Path(argv[1])
