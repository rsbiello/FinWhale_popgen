#!/bin/bash

#------------------------
# Environment info
#------------------------
echo "Running on:"
pwd
hostname
date

#------------------------
# User-configurable variables
#------------------------
# Input directory containing Beagle files (replace with your path)
INPUT_DIR="/path/to/input"  

# Output directory for PCAngsd results (replace with your path)
OUTPUT_DIR="/path/to/output"  

# Docker image name (should be available locally or pulled)
DOCKER_IMAGE="pcangsd:latest"  

# Beagle input file (must exist inside INPUT_DIR)
BEAGLE_FILE="your_dataset.beagle.gz"  

#------------------------
# Run PCAngsd in Docker
#------------------------
docker run \
    -v "${INPUT_DIR}:/vol/data/" \
    -v "${OUTPUT_DIR}:/vol/output/" \
    --user $(id -u):$(id -g) \
    --rm \
    "${DOCKER_IMAGE}" \
    pcangsd -b "/vol/data/${BEAGLE_FILE}" -o "/vol/output/PCAngsd_covmat"

echo "Job finished at:"
date
