#!/bin/bash

#------------------------
# Job info
#------------------------
echo "Job started on $(hostname) at $(date)"
pwd

#------------------------
# Load ANGSD module
#------------------------
module load angsd-0.932

#------------------------
# User-configurable variables
#------------------------
# Number of populations to model
K=2  

# Input Beagle file (from ANGSD)
GL="/path/to/input_beagle_file.beagle.gz"  

# Path to NGSadmix binary (misc folder of ANGSD)
NGSADMIX_BIN="/path/to/angsd/misc/NGSadmix"  

# Working directory for outputs
WD="/path/to/output_directory"  

# Number of threads to use
CPUS=16  

# List of random seeds for multiple runs
SEEDS=(...)

#------------------------
# Prepare output directory
#------------------------
mkdir -p "$WD"
cd "$WD" || exit

#------------------------
# Run NGSadmix for each seed
#------------------------
for seed in "${SEEDS[@]}"; do
    output_prefix="project_K${K}_seed${seed}"
    echo "Running NGSadmix with K=$K and seed $seed..."
    "$NGSADMIX_BIN" -likes "$GL" -K "$K" -P "$CPUS" -seed "$seed" -o "$output_prefix"
done

echo "Job completed at $(date)"
