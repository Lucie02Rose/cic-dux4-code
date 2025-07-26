#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J methylation-giab-model
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-methylation.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-methylation.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate methylation

# Directories
reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-methylation"
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"

# Create output directory
mkdir -p "$output_dir"

echo "Running methylation on all BAM files in $input_dir..."

# Loop through all BAM files in the input directory and subdirectories
for input_bam in "$input_dir"/*/*.bam; do
    # Extract sample name from BAM filename
    sample_name=$(basename "$input_bam" .bam)

    echo "Processing: $input_bam"

    # Run methylation analysis
    aligned_bam_to_cpg_scores --bam "$input_bam" --output-prefix  "$output_dir/${sample_name}-methylation" --threads 16

    echo "Completed: $input_bam"
done

echo "All methylation analyses completed."


