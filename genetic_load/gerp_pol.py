#!/usr/bin/env python3
"""
vcf_gerp_ancestral_counter.py

Count heterozygous and homozygous derived sites for individuals using:
- a VCF file,
- a population list,
- a BED file with ancestral alleles,
- a BED file with GERP scores.

Usage:
    python vcf_gerp_ancestral_counter.py \
        -l population_list.txt \
        -v input.vcf.gz \
        -b ancestral.bed \
        -g gerp.bed \
        -t 2.0 \
        -o output.txt

Author: Your Name
Date: 2025-08-26
"""

import argparse
import re
import gzip

# ------------------------
# Argument parsing
# ------------------------
parser = argparse.ArgumentParser(description="Count derived allele categories per individual")
parser.add_argument('-l', '--populationlist', required=True, type=str,
                    help='File with individual IDs and populations (tab-separated)')
parser.add_argument('-v', '--variationfile', required=True, type=str,
                    help='VCF file with genotypes (can be .gz)')
parser.add_argument('-b', '--ancestralbed', required=True, type=str,
                    help='BED file with ancestral alleles (chrom start end allele)')
parser.add_argument('-g', '--gerpbed', required=True, type=str,
                    help='BED file with GERP scores (score in 5th column)')
parser.add_argument('-t', '--thres', required=True, type=float,
                    help='GERP threshold for deleterious sites')
parser.add_argument('-o', '--output', required=True, type=str,
                    help='Output file name')
args = parser.parse_args()

# ------------------------
# Load populations
# ------------------------
indDict = {}
popList = []

with open(args.populationlist) as f:
    for line in f:
        ind, pop = line.rstrip().split("\t")
        indDict[ind] = {
            'p': pop,
            'calls': 0,
            'heterozygous_sites': 0,
            'homozygous_der_sites': 0
        }
        if pop not in popList:
            popList.append(pop)

popList.sort()
print(f"Loaded {len(indDict)} individuals from {len(popList)} populations: {popList}")

# ------------------------
# Load ancestral alleles
# ------------------------
ancDict = {}
with open(args.ancestralbed) as f:
    for line in f:
        chrom, start, end, allele = line.rstrip().split("\t")
        if chrom not in ancDict:
            ancDict[chrom] = {}
        ancDict[chrom][end] = {'anc': allele}

# ------------------------
# Load GERP scores above threshold
# ------------------------
with open(args.gerpbed) as f:
    for line in f:
        tabs = line.rstrip().split("\t")
        if len(tabs) < 4:
            continue
        chrom, start, end, score = tabs[:4]
        score = float(score)
        if score >= args.thres and chrom in ancDict and end in ancDict[chrom]:
            ancDict[chrom][end]['gerp'] = score

# ------------------------
# Open VCF and process sites
# ------------------------
vcf_open = gzip.open if args.variationfile.endswith('.gz') else open

c_all = c_ok = c_miss = 0
header = []

with vcf_open(args.variationfile, 'rt') as vcf:
    for line in vcf:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
            continue

        tabs = line.rstrip().split("\t")
        c_all += 1
        chrom, pos = tabs[0], tabs[1]

        if chrom in ancDict and pos in ancDict[chrom] and 'gerp' in ancDict[chrom][pos]:
            anc = ancDict[chrom][pos]['anc']
            alt = tabs[4]
            derived_allele = alt if anc == tabs[3] else tabs[3] if anc == alt else None
            if not derived_allele:
                continue

            for i in range(9, len(header)):
                gt = tabs[i].split(":")[0]
                if gt in ("0/1", "1/0"):
                    indDict[header[i]]['heterozygous_sites'] += 1
                elif gt == "1/1" and derived_allele == tabs[4]:
                    indDict[header[i]]['homozygous_der_sites'] += 1
                indDict[header[i]]['calls'] += 1
            c_ok += 1
        else:
            c_miss += 1

# ------------------------
# Write output
# ------------------------
with open(args.output, 'w') as out:
    out.write("IND\tPOP\tHETEROZYGOUS_SITES\tHOMOZYGOUS_DERIVED_SITES\tSUM\tWEIGHTED_SUM\n")
    for pop in popList:
        for ind in sorted(indDict.keys()):
            if indDict[ind]['p'] == pop:
                het = indDict[ind]['heterozygous_sites']
                hom_der = indDict[ind]['homozygous_der_sites']
                total = het + hom_der
                weighted = het + 2*hom_der
                out.write(f"{ind}\t{pop}\t{het}\t{hom_der}\t{total}\t{weighted}\n")

print(f"Processed {c_all} sites; {c_ok} written, {c_miss} skipped due to missing ancestral/GERP")
