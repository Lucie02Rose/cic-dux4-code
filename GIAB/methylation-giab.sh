#!/bin/bash
#BSUB -n 16
#BSUB -M 150000
#BSUB -R 'span[hosts=1] select[mem>150000] rusage[mem=150000]'
#BSUB -q long
#BSUB -J methylation-giab-model
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate methylation

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
output_dir="/lustre/scratch126/casm/team274sb/lr26/methylation_giab"
input_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio"

# Create output directory
mkdir -p "$output_dir"

echo "Running methylation on all BAM files in $input_dir..."

# Loop through all BAM files in the input directory and subdirectories
for input_bam in "$input_dir"/*.bam; do
    # Extract sample name from BAM filename
    sample_name=$(basename "$input_bam" .bam)

    echo "Processing: $input_bam"

    # Run methylation analysis
    aligned_bam_to_cpg_scores --bam "$input_bam" --output-prefix  "$output_dir/${sample_name}-methylation" --threads 16

    echo "Completed: $input_bam"
done

echo "All methylation analyses completed."


