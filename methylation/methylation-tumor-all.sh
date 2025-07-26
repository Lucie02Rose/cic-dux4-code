#!/bin/bash
#BSUB -n 16
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q basement
#BSUB -J methylation-tumor-all-countmode
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate methylation

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
output_dir="/lustre/scratch126/casm/team274sb/lr26/methylation_merged_tumor"
input_bam="/lustre/scratch126/casm/team274sb/lr26/pbmm2-tumor-all-from-bams/tumor_final_merged_wtags.bam"

# Create output directory
mkdir -p "$output_dir"

echo "Running methylation on all BAM files..."

# Loop through all BAM files in subdirectories
sample_name=$(basename "$input_bam" .bam)  # Extract sample name

echo "Processing: $input_bam"

aligned_bam_to_cpg_scores --bam "$input_bam" --output-prefix "$output_dir/${sample_name}-methylation" --threads 16 --pileup-mode count

echo "Completed: $input_bam"


echo "All methylation analyses completed."

