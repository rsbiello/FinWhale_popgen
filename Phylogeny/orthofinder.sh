#!/bin/bash -e

# Print basic job info
pwd; hostname; date

# ================================
# Activate Conda Environment
# ================================
conda activate env

# ================================
# Define Input Variables
# ================================

# Working directory containing FASTA files
# Each species should have one FASTA file with protein sequences
WD=/path/to/orthofinder/input

# Move to working directory
cd $WD

# ================================
# Run OrthoFinder
# ================================
# -t : threads for sequence search
# -a : threads for tree inference
# -M : multiple sequence alignment method
# -S : sequence search program (DIAMOND in ultra-sensitive mode)
# -T : tree inference program (RAxML-NG)
# -f : input folder with protein FASTA files
orthofinder \
  -t 8 \
  -a 8 \
  -M msa \
  -S diamond_ultra_sens \
  -T raxml-ng \
  -f $WD

