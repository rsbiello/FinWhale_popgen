#!/bin/bash

# Basic job info
pwd; hostname; date

# Load Chromap (make sure module system is configured on your HPC cluster)
module load chromap-0.2.4

# ================================
# Define Input/Output Variables
# ================================

# Reference genome (FASTA file)
REF=/path/to/reference/genome.fna

# Input Hi-C sequencing data (FASTQ files)
R1=/path/to/input/sample_R1.fastq.gz
R2=/path/to/input/sample_R2.fastq.gz

# Prefix for output files
ID=sample_project_name

# Working directory for Chromap output
WD=/path/to/output/directory

# Navigate to working directory
cd $WD

# ================================
# Chromap Workflow
# ================================

# (Optional) Build Chromap index if not already built
# This creates an index based on the reference genome
# Uncomment and run once before alignment
# chromap -i -r $REF -o $ID

# Run Chromap with Hi-C preset and generate outputs in multiple formats

# Generate Hi-C pairs format
chromap --preset hic -x $ID -r $REF -1 $R1 -2 $R2 \
  --remove-pcr-duplicates --pairs -o $ID.pairs

# Generate BED format output
chromap --preset hic -x $ID -r $REF -1 $R1 -2 $R2 \
  --remove-pcr-duplicates --BED -o $ID.bed

# Generate SAM format output
chromap --preset hic -x $ID -r $REF -1 $R1 -2 $R2 \
  --remove-pcr-duplicates --SAM -o $ID.sam

