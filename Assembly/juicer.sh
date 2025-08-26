#!/bin/bash

# Print basic job info
pwd; hostname; date

# ================================
# Load required software modules
# ================================
module load yahs-1.2
module load juicertools-2.16
module load samtools-1.11

# ================================
# Define Input/Output Variables
# ================================

# YAHS output files
BIN=/path/to/yahs/output/yahs.out.bin
AGP=/path/to/yahs/output/yahs.out_scaffolds_final.agp

# Reference genome
REF=/path/to/reference/genome.fna

# Working directory for Juicer output
WD=/path/to/juicer/output

# ================================
# Prepare Output Directory
# ================================
mkdir -p $WD
cd $WD

# ================================
# Step 1: Index Reference Genome
# ================================
# This creates a FASTA index (.fai) required by Juicer
samtools faidx $REF

# ================================
# Step 2: Run Juicer Pre-processing
# ================================
# The "juicer pre" step prepares input files for assembly visualization
juicer pre -a -o out_JBAT $BIN $AGP $REF.fai > out_JBAT.log 2>&1

# ================================
# Step 3: Generate .hic File
# ================================
# Convert pre-processed results into a .hic file (for 3D genome visualization in Juicebox)
(java -jar -Xmx250G /path/to/juicertools/juicer_tools_2.16.00.jar pre \
  out_JBAT.txt out_JBAT.hic.part \
  <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) \
  && (mv out_JBAT.hic.part out_JBAT.hic)
