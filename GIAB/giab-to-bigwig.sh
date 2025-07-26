#!/bin/bash
#BSUB -n 32
#BSUB -M 80000
#BSUB -R 'span[hosts=1] select[mem>80000] rusage[mem=80000]'
#BSUB -q basement
#BSUB -J pbmm2-giab
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate deeptools

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
input_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Loop through all BAM files in the input directory
for input_bam in "$input_dir"/*.bam; do
    # Extract base name of the BAM file (without extension)
    base_name=$(basename "$input_bam" .bam)
    
    # Define output BAM file path
    output_bw="$output_dir/${base_name}_pbmm2.bw"

    # Print which file is being processed
    echo "BAM to BW"

    # Run the pbmm2 alignment
    bamCoverage -b "$input_bam" -o "$output_bw" --binSize 50 --normalizeUsing RPKM --effectiveGenomeSize 3100000000 --extendReads --ignoreDuplicates

    # Print completion message
    echo "Conversion completed: $output_bw"
done
