#!/bin/bash
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J cnv
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-hificnv.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-hificnv.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate hificnv || { echo "Failed to activate Conda environment"; exit 1; }

reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
bam_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/blood_1C01_hifi_reads_pbmm2.bam"
output_vcf_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-hificnv-tumor-all"

cd "$output_vcf_dir"

hificnv \
    --bam "$bam_file" \
    --ref "$reference" \
    --threads 16 \
    --output-prefix blood_cnv
