#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J 1A01
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-1A01.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-1A01.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_bam="/lustre/scratch126/cellgen/behjati/lr26/PacBio/tumor_1A01_hifi_reads.bam"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"

# Define output BAM file name
base_name=$(basename "$input_bam" .bam)
output_bam="$output_dir/${base_name}_pbmm2.bam"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Run pbmm2 Alignment
echo "Aligning $input_fastq to $reference..."
pbmm2 align "$reference" "$input_bam" "$output_bam" --preset HIFI --sort -j 16 --unmapped

echo "pbmm2 alignment completed: $output_bam"
