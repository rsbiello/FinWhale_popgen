#!/usr/bin/env python3
"""
polarize_snps.py

Assign the most likely ancestral state for each site based on outgroup genotypes.

Usage:
    python polarize_snps.py -i input_table.txt -o output_prefix [-c 5]
"""

import argparse
import re
import sys

#------------------------
# Command-line arguments
#------------------------
parser = argparse.ArgumentParser(
    description="Assign the most likely ancestral state from a table of sites and genotypes."
)
parser.add_argument(
    '-i', '--infile', 
    help='Table with sites and genotypes', 
    type=str, required=True
)
parser.add_argument(
    '-c', '--colStart', 
    help='Column where genotypes start (1-based; default: 5)', 
    type=int, default=5
)
parser.add_argument(
    '-o', '--oprefix', 
    help='Output prefix', 
    type=str, required=True
)
args = parser.parse_args()

col_start = args.colStart - 1  # convert to 0-based indexing

#------------------------
# Main processing
#------------------------
print(f"Processing file {args.infile}; Outgroups starting at column {args.colStart} (1-based)")

output_file = f"{args.oprefix}.ancestral.txt"

with open(output_file, 'w') as out, open(args.infile, 'r') as infile:
    for line in infile:
        if line.startswith("#"):
            continue  # skip header or comment lines
        line = line.rstrip()
        tabs = line.split("\t")

        # Initialize nucleotide counts
        nt_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        total = 0.0

        # Count alleles from outgroup columns
        for col in tabs[col_start:]:
            if col != "N":
                total += 1
                alleles = col.split("/")
                if len(alleles) > 1:
                    nt_counts[alleles[0]] += 0.5
                    nt_counts[alleles[1]] += 0.5
                else:
                    nt_counts[alleles[0]] += 1

        # Determine ancestral allele
        ancestral = "N"
        support = 0.0
        sorted_nt = sorted(nt_counts, key=nt_counts.get, reverse=True)
        best, second = sorted_nt[:2]

        if nt_counts[best] == nt_counts[second]:
            sys.stderr.write(
                f"{tabs[0]}:{tabs[1]}, equal support for {best} and {second}\n"
            )
        else:
            ancestral = best
            support = nt_counts[best] / total if total > 0 else 0

        # Write output
        out.write(
            "\t".join(tabs[:col_start]) + f"\t{ancestral}\t{support:.2f}\t{int(total)}\n"
        )
