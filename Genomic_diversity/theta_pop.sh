#!/bin/bash
#================================================================
# ANGSD Theta Pipeline: CDS and non-CDS regions
# Author: Roberto Biello
# Date: 2025
# Description: Extract CDS coordinates, create callable regions,
# compute SAFs and estimate theta for SOC population (can extend to others)
#================================================================

echo "Job started on $(hostname) at $(date)"
pwd

#------------------------
# Load modules
#------------------------
module load bedtools-2.30
module load angsd-0.932

#------------------------
# User paths and variables
#------------------------
GFF3="/path/to/file.gff3"
CALLABLE_BED="/path/to/callable.bed"
REF="/path/to/ref.fa"
CHR_LIST="/path/to/chr_list"
BAM_LIST="/path/to/bam_list"
OUT_DIR="/path/to/output/"
ANGSD_MISC="/software/ngs/angsd-0.932/misc"
ID="POP1"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

#------------------------
# Step 1: Extract CDS coordinates
#------------------------
echo "Extracting CDS coordinates..."
cat "$GFF3" | grep "CDS" | awk -F'\t' '{print $1, $4-1, $5}' \
  | sort -k1,1 -k2,2n \
  | sed 's/ /\t/g' \
  > cds_coordinates_sorted_tab.bed

# Merge overlapping CDS regions
bedtools merge -i cds_coordinates_sorted_tab.bed > cds_coordinates_merged.bed

# Intersect with callable regions to get CDS-only callable
bedtools intersect -a cds_coordinates_merged.bed -b "$CALLABLE_BED" > callable_cds.bed

# Count callable bases
awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' callable_cds.bed

# Prepare ANGSD sites file
awk '{print $1"\t"$2+1"\t"$3}' callable_cds.bed > callable_cds.angsd.file
angsd sites index callable_cds.angsd.file

#------------------------
# Step 2: Callable regions without CDS
#------------------------
echo "Generating callable regions excluding CDS..."
bedtools subtract -a "$CALLABLE_BED" -b cds_coordinates_merged.bed > callable_no_cds.bed

# Count bases
awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' callable_no_cds.bed

awk '{print $1"\t"$2+1"\t"$3}' callable_no_cds.bed > callable_no_cds.angsd.file
angsd sites index callable_no_cds.angsd.file

# Clean temporary files
rm cds_coordinates*.bed

#------------------------
# Step 3: Compute SAF using ANGSD
#------------------------
echo "Computing SAF for CDS and non-CDS regions..."
angsd -bam "$BAM_LIST" -out saf.${ID}.cds \
      -minMapQ 20 -sites callable_cds.angsd.file -ref "$REF" -rf "$CHR_LIST" -minQ 20 \
      -doSaf 1 -GL 2 -P 16 -anc "$REF"

angsd -bam "$BAM_LIST" -out saf.${ID}.no_cds \
      -minMapQ 20 -sites callable_no_cds.angsd.file -ref "$REF" -rf "$CHR_LIST" -minQ 20 \
      -doSaf 1 -GL 2 -P 16 -anc "$REF"

#------------------------
# Step 4: Estimate Theta using realSFS and thetaStat
#------------------------
echo "Estimating theta..."
for region in cds no_cds; do
  SAF_FILE="saf.${ID}.${region}.saf.idx"
  OUT_FILE="baphy_${ID}.${region}"

  # Compute folded SFS
  "$ANGSD_MISC/realSFS" "$SAF_FILE" -P 16 -fold 1 > "$OUT_FILE.sfs"

  # Convert SAF to theta
  "$ANGSD_MISC/realSFS" saf2theta "$SAF_FILE" -outname "$OUT_FILE" -sfs "$OUT_FILE.sfs" -fold 1

  # Compute theta statistics per scaffold
  "$ANGSD_MISC/thetaStat" do_stat "$OUT_FILE.thetas.idx"
done

echo "Pipeline completed at $(date)"
