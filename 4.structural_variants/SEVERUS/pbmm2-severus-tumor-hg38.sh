#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J severus-hg38-tum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-severus-hg38.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-severus-hg38.err

### activate the severus conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate severus_env

### define the input, output and reference ###
reference="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.fa"
input_bam_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor_all_4_hifi_reads"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-hg38"
### create the output directory ###
mkdir -p "$output_dir"
### same settings as with the previous somatic callers ###
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
        --vaf-thr 0.01 \
        --min-mapq 10 \
        --min-sv-size 20 \
        --bp-cluster-size 100 \
        --output-read-ids \
        --single-bp

    echo "Severus SV calling completed for $base_name. Results in $severus_sample_dir"
done
