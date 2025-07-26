#!/bin/bash
#BSUB -n 32
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J pbmm2-nanomonsv-match
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-nanomon.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-nanomon.err

set -euo pipefail
echo "Job started on $(date)"
echo "Running on $(hostname)"

# Step 0: Set working directory
cd /lustre/scratch126/cellgen/behjati/lr26 || exit 1
echo "Working directory: $(pwd)"

# Load conda
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate nanomonsv

# Define paths
PARSE_DIR="$(pwd)/nanomonsv-parse_1"
OUTPUT_DIR="$(pwd)/nanomonsv-results-matched_1"
REFERENCE="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
SAMPLE_FASTQ="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_4_hifi_reads.fastq.gz"
CONTROL_FASTQ="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/blood_1C01_hifi_reads.fastq.gz"
SAMPLE_NAME="tumor-all"
CONTROL_NAME="blood"

# Create directories
mkdir -p fastq/ bam/$SAMPLE_NAME/ bam/$CONTROL_NAME/ output/$SAMPLE_NAME/ output/$CONTROL_NAME/ "$OUTPUT_DIR" "$PARSE_DIR"

# Step 1: Decompress FASTQ
echo "Decompressing FASTQ files..."
gunzip -c "$SAMPLE_FASTQ" > "fastq/$SAMPLE_NAME.fastq"
gunzip -c "$CONTROL_FASTQ" > "fastq/$CONTROL_NAME.fastq"

# Step 2: Align reads
echo "Aligning $SAMPLE_NAME..."
minimap2 -ax map-hifi -t 32 "$REFERENCE" "fastq/$SAMPLE_NAME.fastq" | samtools view -Shb - > "bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam"

echo "Sorting and indexing $SAMPLE_NAME..."
samtools sort -@ 32 -m 2G "bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam" -o "bam/$SAMPLE_NAME/$SAMPLE_NAME.bam"
samtools index "bam/$SAMPLE_NAME/$SAMPLE_NAME.bam"

echo "Aligning $CONTROL_NAME..."
minimap2 -ax map-hifi -t 32 "$REFERENCE" "fastq/$CONTROL_NAME.fastq" | samtools view -Shb - > "bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam"

echo "Sorting and indexing $CONTROL_NAME..."
samtools sort -@ 32 -m 2G "bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam" -o "bam/$CONTROL_NAME/$CONTROL_NAME.bam"
samtools index "bam/$CONTROL_NAME/$CONTROL_NAME.bam"

# Step 3: Parse SVs
echo "Parsing SVs for $SAMPLE_NAME..."
nanomonsv parse "bam/$SAMPLE_NAME/$SAMPLE_NAME.bam" "output/$SAMPLE_NAME/$SAMPLE_NAME"

echo "Parsing SVs for $CONTROL_NAME..."
nanomonsv parse "bam/$CONTROL_NAME/$CONTROL_NAME.bam" "output/$CONTROL_NAME/$CONTROL_NAME"

# Step 4: Get SVs with control sample
echo "Running SV detection with NanoMonsv..."
nanomonsv get \
    "output/$SAMPLE_NAME/$SAMPLE_NAME" \
    "bam/$SAMPLE_NAME/$SAMPLE_NAME.bam" \
    "$REFERENCE" \
    --control_prefix "output/$CONTROL_NAME/$CONTROL_NAME" \
    --control_bam "bam/$CONTROL_NAME/$CONTROL_NAME.bam" \
    --single_bnd \
    --use_racon \
    --processes 32 \
    --min_tumor_variant_read_num 2 \
    --min_tumor_VAF 0.01 \
    --max_control_variant_read_num 3 \
    --cluster_margin_size 100 \
    --median_mapQ_thres 25 \
    --max_overhang_size_thres 30 \
    --var_read_min_mapq 10 \
    > "$OUTPUT_DIR/nanomonsv_get.log" 2>&1

# Step 5: Cleanup
echo "Cleaning up temporary files..."
rm "fastq/$SAMPLE_NAME.fastq" "fastq/$CONTROL_NAME.fastq"
rm "bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam" "bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam"

echo "Job completed on $(date)"
