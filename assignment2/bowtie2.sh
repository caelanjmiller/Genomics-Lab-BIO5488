#!/bin/bash

INDEX_PREFIX="$1"
READS="$2"
SAM_OUTPUT="$3"
# If MAC OS
THREADS=$(sysctl -n hw.ncpu)
# All else Linux
# THREADS=$(nproc --all)

bowtie2 --threads "$THREADS" -x "$INDEX_PREFIX" -U "$READS" -S "$SAM_OUTPUT" 2> alignment_report.txt 
