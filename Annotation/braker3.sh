#!/bin/bash -e

# Logging
pwd; hostname; date

# === Load environment ===
conda activate braker3

# Ensure GeneMark is in PATH (adjust to your environment)
export PATH="/path/to/genemark/bin:$PATH"

# === Input files ===
GENOME=/path/to/genome/genome.fa.masked
PROT_VER=/path/to/proteins/vertebrata_proteins.fa
PROT_BAL=/path/to/proteins/species_proteins.fa
BAM=/path/to/rna_alignments/sample_merged.junc_filt.bam

# === Working directory & species ID ===
ID=example_species.braker3
WD=/path/to/output/braker3_output
mkdir -p $WD

# === Run BRAKER3 ===
braker.pl \
  --genome=$GENOME \
  --workingdir=$WD \
  --PROTHINT_PATH=/path/to/prothint/bin \
  --threads 16 \
  --gff3 \
  --bam=$BAM \
  --prot_seq=$PROT_VER,$PROT_BAL \
  --species=$ID

date

