#!/bin/bash
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J sawfish
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sawfish.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sawfish.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate sawfish || { echo "Failed to activate Conda environment"; exit 1; }

reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.fa"
bam_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/1C01/blood_1C01_hifi_reads_pbmm2.bam"
output_vcf_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-blood"

# Extract base name from BAM file
base_name=$(basename "$bam_file" .bam)

# Define unique output directories
discover_dir="$output_vcf_dir/${base_name}_discover"
joint_call_dir="$output_vcf_dir/${base_name}_joint_call"

# Remove existing directories if they exist
[ -d "$discover_dir" ] && rm -rf "$discover_dir"
[ -d "$joint_call_dir" ] && rm -rf "$joint_call_dir"

echo "Running Sawfish discover on $bam_file using $reference..."

# Step 1: Discover SVs
sawfish discover --bam "$bam_file" --output-dir "$discover_dir" --ref "$reference" --threads 16

echo "Running Sawfish joint-call for $base_name..."

# Step 2: Joint calling
sawfish joint-call --sample "$discover_dir" --output-dir "$joint_call_dir" --threads 16

echo "Sawfish processing completed for: $base_name"
