#!/bin/bash
#BSUB -n 40
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J t2t-alignment-single
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/bamscale.%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/bamscale.%J.e

set -euo pipefail

echo "Activating environment"
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bamscale

# Set parameters
threads=40
roi_bed="/lustre/scratch126/casm/team274sb/lr26/population_dux/dux.bed"
output_dir="/lustre/scratch126/casm/team274sb/lr26/population_dux"
bam_dir="$output_dir"
bam_out="$output_dir/bamscale_output"

mkdir -p "$bam_out"

# Run BAMscale in parallel (1 thread per BAM, 40 BAMs at a time)
ls "$bam_dir"/*.bam | parallel -j $threads '
  sample=$(basename {} .bam);
  echo "Processing $sample";
  BAMscale cov --bed "'$roi_bed'" \
    --bam {} \
    --outdir "'$bam_out'" \
    --threads 1 \
    --mapq 10 \
    --frag \
    --prefix "dux4_${sample}"
'
