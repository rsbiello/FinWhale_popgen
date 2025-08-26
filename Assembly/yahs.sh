#!/bin/bash

# Basic job info
pwd; hostname; date

# Load YAHS (make sure module system is configured on your HPC cluster)
module load yahs-1.2

# ================================
# Define Input/Output Variables
# ================================

# Reference genome (FASTA file)
REF=/path/to/reference/genome.fna

# Input BAM file (aligned Hi-C reads, typically from Chromap or BWA-MEM)
BAM=/path/to/input/sample.bam

# ================================
# Run YAHS
# ================================

# The -e flag specifies the restriction enzyme recognition site used for Hi-C data
# (example: "GATC" for DpnII; adjust based on your experiment)
yahs -e "GATC" $REF $BAM
