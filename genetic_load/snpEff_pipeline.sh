#!/bin/bash

################################################################################
# SNP EFFECT PIPELINE
# This script performs:
# 1) Extract variants from VCF
# 2) Run SnpEff annotation
# 3) Extract coding regions
# 4) Separate variants into MODERATE, LOW, HIGH
################################################################################

# === CONFIGURATION ===
VCF_RAW="/path/to/file.vcf.gz"
SNPEFF_HOME="/path/to/snpEff"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
DATA_DIR="${SNPEFF_HOME}/data"
GENOME="species"

WORKDIR="$(pwd)/output/"
OUT_PREFIX="species_snpeff"

mkdir -p $WORKDIR
cd $WORKDIR

################################################################################
# STEP 1: Extract non-missing variants
echo "STEP 1: Extracting variants..."
zcat $VCF_RAW | awk '$5!="."' | gzip -c > ${OUT_PREFIX}.variants.vcf.gz
echo "Variants extracted: ${OUT_PREFIX}.variants.vcf.gz"

################################################################################
# STEP 2: Run SnpEff annotation
echo "STEP 2: Running SnpEff..."
# Uncomment the next line to build the genome database (only needed once)
# java -Xmx4g -jar ${SNPEFF_HOME}/snpEff.jar build -gff3 -c $SNPEFF_CONFIG -nodownload -dataDir $DATA_DIR -noCheckCds -noCheckProtein -v $GENOME

java -Xmx4g -jar ${SNPEFF_HOME}/snpEff.jar $GENOME ${OUT_PREFIX}.variants.vcf.gz | gzip -c > ${OUT_PREFIX}.snpeff.vcf.gz
echo "SnpEff annotation completed: ${OUT_PREFIX}.snpeff.vcf.gz"

################################################################################
# STEP 3: Extract coding regions from annotation GTF
echo "STEP 3: Extract CDS regions..."
GTF="${DATA_DIR}/${GENOME}/genes.gtf.gz"
zcat $GTF | awk '$3=="CDS"' > ${OUT_PREFIX}.cds.gtf
bedtools intersect -header -a ${OUT_PREFIX}.snpeff.vcf.gz -b ${OUT_PREFIX}.cds.gtf -u > ${OUT_PREFIX}.snpeff.cds.vcf
echo "CDS VCF created: ${OUT_PREFIX}.snpeff.cds.vcf"

################################################################################
# STEP 4: Separate variants by impact
echo "STEP 4: Splitting by variant impact..."
bcftools view -i 'INFO/ANN[*] ~ "MODERATE"' ${OUT_PREFIX}.snpeff.cds.vcf > ${OUT_PREFIX}.snpeff.cds.MODERATE.vcf
bcftools view -i 'INFO/ANN[*] ~ "LOW"'      ${OUT_PREFIX}.snpeff.cds.vcf > ${OUT_PREFIX}.snpeff.cds.LOW.vcf
bcftools view -i 'INFO/ANN[*] ~ "HIGH"'     ${OUT_PREFIX}.snpeff.cds.vcf > ${OUT_PREFIX}.snpeff.cds.HIGH.vcf
echo "Variants split into MODERATE, LOW, HIGH."

echo "Pipeline finished!"
