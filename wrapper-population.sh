#!/bin/bash
#BSUB -n 32
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q long
#BSUB -J wgs-t2t-alignment-single
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/pop.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/pop.e

set -euo pipefail

if [[ -z "${SAMPLE_BAM:-}" ]]; then
    echo "Please set SAMPLE_BAM=/full/path/to/sample.bam"
    exit 1
fi

echo "Running population-dux.sh on $SAMPLE_BAM"
bash population-dux.sh "$SAMPLE_BAM"
