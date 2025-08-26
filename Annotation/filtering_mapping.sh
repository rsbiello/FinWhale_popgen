#!/bin/bash

# Print basic job info
pwd; hostname; date

# ================================
# Load required software
# ================================
# Assumes fastp is installed in $SD, and Portcullis is available in the conda env
conda activate portcullis   # Load conda environment with Portcullis
module load hisat2-2.2.1
module load samtools-1.11

# ================================
# Define Input/Output Variables
# ================================

SD=/path/to/software                     # Directory containing fastp
REF=/path/to/reference/genome.fna        # Reference genome (FASTA)
WD=/path/to/working/directory            # Working directory for RNA-seq analysis
REF_PREFIX=$WD/genome_index              # Prefix for Hisat2 index
LIST=$WD/sample_IDs.txt                  # File containing list of sample IDs (one per line)
FQ=/path/to/raw/fastq                    # Directory containing raw FASTQ files
THREADS=8

OUT_TRIM=$WD/trimmed                     # Output: trimmed reads
OUT_BAM=$WD/BAMs                         # Output: BAM files
OUT_PORT=$WD/portcullis                  # Output: Portcullis results

# Create output directories
mkdir -p $OUT_TRIM $OUT_BAM $OUT_PORT

# ================================
# Step 1: Quality Filtering (fastp)
# ================================
echo "Step 1: Quality filtering with fastp..."
for ID in $(cat $LIST); do
    R1=$FQ/${ID}_1.fastq.gz
    R2=$FQ/${ID}_2.fastq.gz
    OUT1=$OUT_TRIM/${ID}_R1.filt.fq.gz
    OUT2=$OUT_TRIM/${ID}_R2.filt.fq.gz

    echo "Processing $ID..."
    $SD/fastp -q 20 -l 75 -i $R1 -I $R2 -o $OUT1 -O $OUT2 || { echo "fastp failed for $ID" >&2; exit 1; }
done

# ================================
# Step 2: Build Hisat2 Index
# ================================
echo "Step 2: Building Hisat2 index..."
hisat2-build -p $THREADS $REF $REF_PREFIX || { echo "Hisat2 index build failed" >&2; exit 1; }

# ================================
# Step 3: Align Reads and Process BAMs
# ================================
echo "Step 3: Aligning reads with Hisat2..."
for ID in $(cat $LIST); do
    R1=$OUT_TRIM/${ID}_R1.filt.fq.gz
    R2=$OUT_TRIM/${ID}_R2.filt.fq.gz
    SAM=$OUT_BAM/$ID.sam
    BAM_OUT=$OUT_BAM/$ID.bam
    SORTED_BAM=$OUT_BAM/$ID.sorted.bam
    TMP_DIR=$OUT_BAM/$ID.tmp

    echo "Aligning $ID..."
    hisat2 -q --dta-cufflinks -p $THREADS -x $REF_PREFIX -1 $R1 -2 $R2 -S $SAM || { echo "Hisat2 failed for $ID" >&2; exit 1; }

    echo "Converting SAM to BAM..."
    samtools view -b -o $BAM_OUT $SAM || { echo "SAM to BAM failed for $ID" >&2; exit 1; }

    echo "Sorting BAM..."
    mkdir -p $TMP_DIR
    samtools sort -o $SORTED_BAM -T $TMP_DIR/aln.sorted --threads $THREADS -m 10G $BAM_OUT || { echo "Sorting BAM failed for $ID" >&2; exit 1; }

    echo "Indexing BAM..."
    samtools index -c $SORTED_BAM || { echo "Indexing BAM failed for $ID" >&2; exit 1; }
done

# ================================
# Step 4: Index Reference FASTA
# ================================
echo "Step 4: Indexing reference FASTA..."
samtools faidx $REF || { echo "FASTA indexing failed" >&2; exit 1; }

# ================================
# Step 5: Run Portcullis Pipeline
# ================================
echo "Step 5: Running Portcullis..."
for ID in $(cat $LIST); do
    PREP_DIR=$OUT_PORT/1_prep/$ID
    JUNC_DIR=$OUT_PORT/2_junc/$ID
    FILT_DIR=$OUT_PORT/3_filt/$ID
    BAM_FILT=$OUT_PORT/4_BAM_filt/$ID.junc_filt.bam

    echo "Preparing Portcullis for $ID..."
    portcullis prep --use_csi -o $PREP_DIR $REF $OUT_BAM/$ID.sorted.bam || { echo "Portcullis prep failed for $ID" >&2; exit 1; }

    echo "Running junction detection..."
    portcullis junc --use_csi -t $THREADS -o $JUNC_DIR $PREP_DIR || { echo "Junction detection failed for $ID" >&2; exit 1; }

    echo "Filtering junctions..."
    portcullis filt -t $THREADS -o $FILT_DIR $PREP_DIR $JUNC_DIR.junctions.tab || { echo "Filtering failed for $ID" >&2; exit 1; }

    echo "Filtering BAM..."
    portcullis bamfilt --use_csi -o $BAM_FILT $FILT_DIR.pass.junctions.tab $PREP_DIR/portcullis.sorted.alignments.bam || { echo "BAM filtering failed for $ID" >&2; exit 1; }
done

echo "Pipeline completed successfully!"
