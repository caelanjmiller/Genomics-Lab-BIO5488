#!/bin/bash

CHAR_SEARCH=$1
FASTA=$2
# Skip 1st line of FASTA file (header), capture the number of occurrences of provided character & count it
COUNT=$(tail -n +2 "${FASTA}" | grep -o "${CHAR_SEARCH}" | wc -l)
echo "${CHAR_SEARCH}" appears "${COUNT}" times in "${FASTA}"