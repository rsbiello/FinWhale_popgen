#!/bin/bash -e

############################################
# Paths
############################################
REF=/jic/scratch/projects/CA792-D11-B/bal/reference/Baphy.unife.v2.scaffolds.fa
WORK=/jic/scratch/projects/CA792-D11-B/bal/masking

############################################
# 1. RepeatMasker with mammalia db (softmask)
############################################
source package repeatmasker-4.0.7
WD1=$WORK/repeatmasker/1run_mammalia
mkdir -p $WD1 && cd $WD1

RepeatMasker -pa 8 -e ncbi -s -species mammalia -xsmall -gff $REF

############################################
# 2. Extract taxon-specific repeats (mammalia)
############################################
source package bioperl-1.6.923
perl queryRepeatDatabase.pl \
    -species mammalia > repeatmasker.mammalia.fa

############################################
# 3. RepeatModeler 2.0.2a
############################################
conda activate RepeatModeler
source package recon-1.08

WD2=$WORK/repeatmodeller
mkdir -p $WD2 && cd $WD2

BuildDatabase -name "B_physalus" -engine ncbi $REF
RepeatModeler -engine ncbi -pa 16 -trf_prgm /hpc-home/biello/trf409.linux64 -database B_physalus

############################################
# 4. Merge libraries (ReannTE)
############################################
MASKER=$WD1/repeatmasker.mammalia.fa
MODELLER=$WD2/*/consensi.fa.classified
WD3=$WORK/merged
mkdir -p $WD3 && cd $WD3

ReannTE_MergeFasta.pl \
    -a $MASKER -b $MODELLER \
    -RM /repeatmasker/4.0.7/x86_64/bin/ -CheckLow 80

############################################
# 5. Final RepeatMasker run with merged lib + TRF
############################################
WD4=$WORK/repeatmasker/2run_merged
mkdir -p $WD4 && cd $WD4

LIB=$(ls $WD2/*/consensi.fa.classified.nolow.fa)
RepeatMasker -pa 8 -dir . -trf_prgm trf409.linux64 \
    -xsmall -a -no_is -e ncbi -s -gff -lib $LIB $REF

date
echo "Repeat annotation pipeline finished successfully!"
