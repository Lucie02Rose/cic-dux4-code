#!/bin/bash
### LSF job parameters ###
#BSUB -n 32
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q long
#BSUB -J wgs-pop_wrapper
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-pop.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-pop.e

### failure management ###
set -euo pipefail

### sample bams ###
if [[ -z "${SAMPLE_BAM:-}" ]]; then
    echo "Please set SAMPLE_BAM=/full/path/to/sample.bam"
    exit 1
fi

echo "Running population-dux.sh on $SAMPLE_BAM"
### choose the name of script to be run (there are number of the population dux scripts) ###
bash population-dux.sh "$SAMPLE_BAM"
