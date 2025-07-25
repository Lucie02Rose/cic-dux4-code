#!/bin/bash
#BSUB -n 16
#BSUB -M 75000
#BSUB -R 'span[hosts=1] select[mem>75000] rusage[mem=75000]'
#BSUB -q normal
#BSUB -J pbmm2-sev-1B01
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-severus.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-severus.err

# Activate conda environment for Severus
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate severus_env

# Directories
reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.fa"
input_tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
input_blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/blood_1C01_hifi_reads_pbmm2.bam"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-new"

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Process each BAM file in the directory
echo "Running Severus..."

severus \
   --target-bam "$input_tumor" \
   --control-bam "$input_blood" \
   --out-dir "$output_dir" \
   -t 32 \
   --min-support 2 \
   --vaf-thr 0.01 \
   --min-mapq 10 \
   --bp-cluster-size 100 \
   --min-sv-size 20 \
   --single-bp \
   --output-read-ids

echo "Severus SV call completed"
