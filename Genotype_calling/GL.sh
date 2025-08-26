#!/bin/bash

#------------------------
# Environment info
#------------------------
echo "Running on:"
pwd
hostname
date

#------------------------
# Load ANGSD module
#------------------------
module load angsd-0.932

#------------------------
# User-configurable variables
#------------------------
# List of chromosomes or regions
CHR="/path/to/chromosomes_list"  

# Callable sites mask
MASK="/path/to/callable_sites.file"  

# List of BAM files (one per line)
BAMLIST="/path/to/bam_list.txt"  

# Output ID / prefix
ID="project_name.MAF_05"  

# Working directory for output
WD="/path/to/output_directory"  

# Minimum number of individuals with data at a site
MININD=

# Number of CPUs to use
CPUS=16  

#------------------------
# Prepare output directory
#------------------------
mkdir -p "$WD"
cd "$WD" || exit

#------------------------
# Run ANGSD
#------------------------
angsd \
  -bam "$BAMLIST" \
  -out "$ID" \
  -minMapQ 20 \
  -sites "$MASK" \
  -rf "$CHR" \
  -minQ 20 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -SNP_pval 1e-6 \
  -minInd "$MININD" \
  -GL 2 \
  -doGlf 2 \
  -minMaf 0.05 \
  -P "$CPUS"

echo "Job finished at:"
date
