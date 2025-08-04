#!/bin/bash
### parameters for LSF job ###
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J cnv
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-hificnv-blood.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-hificnv-blood.err

### activate the hificnv conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate hificnv

### define file paths ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
bam_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/blood_1C01_hifi_reads_pbmm2.bam"
output_vcf_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-hificnv"

### make and change to output directory ###
mkdir -p "$output_vcf_dir"
cd "$output_vcf_dir"

### run hificnv on blood ###
hificnv \
    --bam "$bam_file" \
    --ref "$reference" \
    --threads 16 \
    --output-prefix blood_cnv
