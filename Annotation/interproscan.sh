#!/bin/bash

# Logging job info
pwd; hostname; date

# === Load InterProScan environment ===
# Replace this with your module/environment setup
source package interproscan-5.52-86.0

# === Input and output files ===
# PROT: protein sequences (FASTA) to annotate
PROT=/path/to/proteins/my_species.ab_initio.proteins.fa

# OUT: base name for InterProScan output
OUT=my_species.interproscan.Pfam

# WD: working directory for InterProScan results
WD=/path/to/output/interproscan
mkdir -p $WD
cd $WD

# === Step 1: Run InterProScan with Pfam application ===
# -appl Pfam      : run only the Pfam database
# -f XML          : output in XML format
# -i PROT         : input protein FASTA
# -b OUT          : base name for output files
# -dp             : include detailed protein matches
# -goterms        : add GO term annotations
# -pa             : add pathway annotations
interproscan.sh -appl Pfam -f XML -i $PROT -b $OUT -dp -goterms -pa

# === Step 2: Convert XML output to GFF3 format ===
# -mode convert   : convert mode
# -f GFF3         : output format GFF3
# -i OUT.xml      : XML input
# -b OUT          : base name for converted output
interproscan.sh -mode convert -f GFF3 -i $OUT.xml -b $OUT

# End logging
date
echo "InterProScan Pfam annotation finished!"
