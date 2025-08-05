#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 40
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J bamscale
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-bamscale.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-bamscale.e

### activate the bamscale conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bamscale

### define bam scale parameters - dux coordinates have been saved in a bed file ###
roi_bed="/lustre/scratch126/casm/team274sb/lr26/population_dux/dux.bed"
output_dir="/lustre/scratch126/casm/team274sb/lr26/population_dux"
bam_dir="$output_dir"
bam_out="$output_dir/bamscale_output"
### make the output directory ###
mkdir -p "$bam_out"

### run bamscale in parallel 1 per core ###
ls "$bam_dir"/*.bam | parallel -j 40 '
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
