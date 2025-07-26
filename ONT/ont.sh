#!/bin/bash
#BSUB -n 44
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J ont
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-ont.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-ont.err

set -euo pipefail
echo "Job started on $(date)"
echo "Running on $(hostname)"

# Step 0: Set working directory
cd /lustre/scratch126/cellgen/behjati/lr26 || exit 1
echo "Working directory: $(pwd)"

# Load conda
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate base

# Define paths
OUTPUT_DIR="/lustre/scratch126/cellgen/behjati/lr26/ONT/results"
REFERENCE="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
SAMPLE_FASTQ="/lustre/scratch126/cellgen/behjati/lr26/ONT/normal/pass/*.fastq.gz"
SAMPLE_COMB="/lustre/scratch126/cellgen/behjati/lr26/ONT/normal/tumor.fastq.gz"
tmp="/lustre/scratch126/cellgen/behjati/lr26/ONT/results/tmp"
SAMPLE_NAME="tumor"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$tmp"

# Step 1: Decompress FASTQ
echo "Combining FASTQ files..."
zcat $SAMPLE_FASTQ | gzip > "$SAMPLE_COMB"

# Step 2: Align reads
echo "Aligning $SAMPLE_NAME..."
minimap2 -ax map-ont -t 44 "$REFERENCE" "$SAMPLE_COMB" | samtools view -Shb - > "$OUTPUT_DIR/$SAMPLE_NAME.unsorted.bam"

echo "Sorting and indexing $SAMPLE_NAME..."
samtools sort -@ 44 -m 2G "$OUTPUT_DIR/$SAMPLE_NAME.unsorted.bam" -o "$OUTPUT_DIR/$SAMPLE_NAME.bam"
samtools index "$OUTPUT_DIR/$SAMPLE_NAME.bam"

echo "Job completed on $(date)"
