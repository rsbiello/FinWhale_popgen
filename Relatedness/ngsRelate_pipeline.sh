#!/bin/bash

##############################################
# Steps:
#   1. Generate Genotype Likelihoods (GLF) with ANGSD
#   2. Extract frequency column for ngsRelate
#   3. Run ngsRelate to estimate relatedness
#
# Requirements:
#   - ANGSD >= 0.932
#   - ngsRelate
#   - SLURM environment
##############################################

echo ">>> Job started on $(hostname) at $(date)"

# --- Load modules ---
module load angsd-0.932

# --- Input files & parameters ---
CHR_LIST="/path/to/21autosomes_list"
MASK="/path/to/BaPhy.callable.angsd.file"
BAMLIST="/path/to/NoSibs_49i.C.NEW.bam.list"
IND_LIST="/path/to/list_inds"                # one individual per line, same order as BAM list
N_IND=49                                     # number of individuals
MIN_IND=37                                   # min number of individuals with data
ID="NoSibs_49i.C.NEW.MAF_05"                 # prefix for outputs

# --- Directories ---
WORKDIR="/path/to/analyses/relate"
GL_DIR="${WORKDIR}/GL"
OUT_DIR="${WORKDIR}/output"

mkdir -p "$GL_DIR" "$OUT_DIR"

##############################################
# Step 1: Genotype Likelihoods
##############################################
cd "$GL_DIR"

echo ">>> Running ANGSD to compute GLF..."
angsd -bam "$BAMLIST" -out "$ID" \
  -minMapQ 20 -minQ 20 \
  -sites "$MASK" -rf "$CHR_LIST" \
  -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3 \
  -minInd "$MIN_IND" \
  -GL 2 -doGlf 3 -minMaf 0.05 -P 16

# record number of sites retained from log
NSITES=$(zcat ${ID}.mafs.gz | wc -l)
NSITES=$((NSITES - 1)) # subtract header
echo ">>> Sites retained after filtering: $NSITES"

##############################################
# Step 2: Extract allele frequencies
##############################################
echo ">>> Extracting allele frequencies..."
zcat "${ID}.mafs.gz" | cut -f5 | sed 1d > freq

##############################################
# Step 3: Run ngsRelate
##############################################
cd "$OUT_DIR"

SD="/path/to/ngsRelate"

echo ">>> Running ngsRelate..."
"$SD/ngsRelate" \
  -g "${GL_DIR}/${ID}.glf.gz" \
  -n "$N_IND" \
  -L "$NSITES" \
  -f "${GL_DIR}/freq" \
  -z "$IND_LIST" \
  -O relate_results

echo ">>> Pipeline finished at $(date)"
