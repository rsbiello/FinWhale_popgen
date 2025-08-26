#!/bin/bash

echo "Job started on $(hostname) at $(date)"
pwd

#------------------------
# Load modules
#------------------------
module load msmc-tools
module load msmc2-2.1.1
module load samtools-1.11
module load bcftools-1.11
module load java-1.8

#------------------------
# User-configurable paths
#------------------------
REF="/path/to/reference.fa"                     # reference genome
BAM_DIR="/path/to/BAMs"                        # BAM files
SCAFFS="/path/to/autosome_list.txt"            # scaffolds
LIST_IDS="/path/to/list_IDs.txt"               # sample IDs
MSMC_OUTPUT="/path/to/MSMC/output"             # main output dir
BED_DIR="$MSMC_OUTPUT/BEDs"
MASK_DIR="$BED_DIR/masks"
INPUT_DIR="$MSMC_OUTPUT/input"

GATK_JAR="/path/to/GenomeAnalysisTK.jar"      # GATK CallableLoci
JAVA_BIN="/path/to/java"

THREADS=8

#------------------------
# Create required directories
#------------------------
mkdir -p "$MSMC_OUTPUT" "$BED_DIR" "$MASK_DIR" "$INPUT_DIR"

#------------------------
# Step 1: Generate VCFs per scaffold and sample
#------------------------
echo "Step 1: Generating per-scaffold VCFs"
for ID in $(cat "$LIST_IDS"); do
    for SCAFF in $(cat "$SCAFFS"); do
        samtools mpileup -q 20 -Q 20 -C 50 -u -r "$SCAFF" -f "$REF" "$BAM_DIR/$ID.bam" \
        | bcftools call -c -V indels \
        | bamCaller.py 15 "$ID.mask.$SCAFF.bed.gz" \
        | gzip -c > "$MSMC_OUTPUT/$ID.$SCAFF.vcf.gz"
    done
done

#------------------------
# Step 2: Calculate depth and generate callable loci
#------------------------
echo "Step 2: Depth calculation and callable loci"
for ID in $(cat "$LIST_IDS"); do
    samtools depth "$BAM_DIR/$ID.bam" -a -o "$BED_DIR/$ID.depth"
    
    meandepth=$(awk '{sum+=$3} END {print sum/NR}' "$BED_DIR/$ID.depth")
    echo "Mean depth for $ID = $meandepth"
    echo "$meandepth" > "$BED_DIR/$ID.meandepth.txt"

    mindepth=$(awk -v md="$meandepth" 'BEGIN {printf "%.0f\n", md/3}')
    maxdepth=$(awk -v md="$meandepth" 'BEGIN {printf "%.0f\n", md*3}')
    
    echo "Generating callable loci for $ID (min=$mindepth, max=$maxdepth)"
    "$JAVA_BIN" -Xmx12g -jar "$GATK_JAR" -T CallableLoci \
        -R "$REF" -I "$BAM_DIR/$ID.bam" \
        --minBaseQuality 20 --minMappingQuality 20 \
        --minDepth "$mindepth" --maxDepth "$maxdepth" \
        -summary "$BED_DIR/$ID.summary_CallableLoci.txt" \
        -o "$BED_DIR/$ID.bed"
done

#------------------------
# Step 3: Create callable masks per scaffold
#------------------------
echo "Step 3: Creating scaffold-specific callable masks"
for ID in $(cat "$LIST_IDS"); do
    cat "$BED_DIR/$ID.bed" | grep CALLABLE | cut -f1-3 > "$MASK_DIR/callable_mask.$ID.allscaffs.bed"
    for SCAFF in $(cat "$SCAFFS"); do
        grep "$SCAFF" "$MASK_DIR/callable_mask.$ID.allscaffs.bed" > "$MASK_DIR/callable_mask.$ID.$SCAFF.bed"
    done
done

#------------------------
# Step 4: Create mappability masks per scaffold
#------------------------
echo "Step 4: Creating mappability masks"
MAPPABILITY_BED="/path/to/mappability/BaPhy.gem.150_004.unique.wig.bed"
for SCAFF in $(cat "$SCAFFS"); do
    grep -w "$SCAFF" "$MAPPABILITY_BED" | cut -f1-3 > "$MASK_DIR/mappability_mask.$SCAFF.bed"
done

#------------------------
# Step 5: Generate MSMC2 input files
#------------------------
echo "Step 5: Generating MSMC2 input files"
for ID in $(cat "$LIST_IDS"); do
    for SCAFF in $(cat "$SCAFFS"); do
        generate_multihetsep.py \
            --mask="$MASK_DIR/callable_mask.$ID.$SCAFF.bed" \
            --mask="$MASK_DIR/mappability_mask.$SCAFF.bed" \
            "$MSMC_OUTPUT/$ID.$SCAFF.vcf.gz" \
            > "$INPUT_DIR/$ID.$SCAFF.msmc.txt"
    done
done

#------------------------
# Step 6: Run MSMC2
#------------------------
echo "Step 6: Running MSMC2"
for ID in $(cat "$LIST_IDS"); do
    msmc2_linux64bit -t "$THREADS" -o "$MSMC_OUTPUT/MSMC/$ID.msmc" "$INPUT_DIR/$ID.*.msmc.txt"
done

echo "Pipeline finished at $(date)"
