#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-hifi-new
#BSUB -G team274
#BSUB -o /nfs/users/nfs_l/lr26/outputs/%J_pbmm2-hifi-new.out
#BSUB -e /nfs/users/nfs_l/lr26/errors/%J_pbmm2-hifi-new.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="$output_dir/tmp"
reference_fasta="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.fa"
reference_index="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"

# Check if the index already exists
if [ ! -f "$reference_index" ]; then
    echo "Index not found at $reference_index. Creating pbmm2 index..."
    pbmm2 index "$reference_fasta" "$reference_index"
    echo "Indexing completed."
else
    echo "Reference index already exists at $reference_index."
fi

# Create temp directory if it doesn't exist
mkdir -p "$output_dir" "$tmp_dir"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Loop through all BAM files in the input directory
for input_bam in "$input_dir"/*.bam; do
    # Extract base name of the BAM file (without extension)
    base_name=$(basename "$input_bam" .bam)
    
    # Define output BAM file path
    output_bam="$output_dir/${base_name}_pbmm2.bam"

    # Print which file is being processed
    echo "Aligning $input_bam to $reference..."

    # Run the pbmm2 alignment
    pbmm2 align "$reference_index" "$input_bam" "$output_bam" --preset HIFI --sort -j 16 --unmapped

    # Print completion message
    echo "pbmm2 alignment completed: $output_bam"
done
