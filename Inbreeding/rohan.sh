#!/bin/bash

# === Logging job info ===
pwd; hostname; date

# === Load ROHan module/environment ===
module load rohan-1.0

# === Input files (placeholders) ===
GENOME=/path/to/reference/genome.fa
BAM=/path/to/sample1.sorted.dedup.bam
MAP=/path/to/callable_regions.map          # Callable sites file from ANGSD or similar
AUTO=/path/to/autosome_list.txt            # List of autosomes/scaffolds
OUTDIR=/path/to/output/roh/sample1

# === Parameters ===
THREADS=8
OUTPREFIX=sample1
ROHMU=5e-04          # Mutation rate per site
TSTV=2.0             # Transition/transversion ratio
SIZE=50000           # Minimum ROH length in bp

# === Prepare output directory ===
mkdir -p $OUTDIR
cd $OUTDIR

# === Run ROHan ===
rohan \
  -t $THREADS \
  -o $OUTPREFIX \
  --rohmu $ROHMU \
  --map $MAP \
  --auto $AUTO \
  --tstv $TSTV \
  --size $SIZE \
  $GENOME $BAM

# === End logging ===
date
echo "ROHan analysis completed for $OUTPREFIX"
