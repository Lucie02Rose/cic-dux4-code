#!/bin/bash
#BSUB -n 16
#BSUB -M 75000
#BSUB -R 'span[hosts=1] select[mem>75000] rusage[mem=75000]'
#BSUB -q normal
#BSUB -J pbmm2-sev-1A01
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Activate conda environment for Severus
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate severus_env

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_bam_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/1_A01"
output_dir="/lustre/scratch126/casm/team274sb/lr26/severus-results1A01"

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Process each BAM file in the directory
for bam_file in "$input_bam_dir"/*.bam; do
    # Extract base name without extension
    base_name=$(basename "$bam_file" .bam)
    
    # Define Severus output directory for this sample
    severus_sample_dir="$output_dir/${base_name}-severus"
    
    # Create output directory for the sample
    mkdir -p "$severus_sample_dir"

    echo "Running Severus for $base_name..."

    severus --target-bam "$bam_file" \
        --out-dir "$severus_sample_dir" \
        -t 16 \
        --min-support 2 \
        --vaf-thr 0.02 \
        --min-mapq 10 \
        --min-sv-size 20 \
        --bp-cluster-size 100 \
        --output-read-ids \
        --single-bp

    echo "Severus SV calling completed for $base_name. Results in $severus_sample_dir"
done
