#!/bin/bash
#BSUB -n 32
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J pbmm2-nanomonsv-match
#BSUB -G team274  # Fix extra space
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-nanomon.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-nanomon.err

# Activate conda environment for Nanomonsv
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate nanomonsv

# Define paths
PARSE_DIR="/lustre/scratch126/cellgen/behjati/lr26/nanomonsv-parse"
OUTPUT_DIR="/lustre/scratch126/cellgen/behjati/lr26/nanomonsv-results-matched"
REFERENCE="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.fa"
SAMPLE_FASTQ="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_hifi_reads.fastq.gz"
CONTROL_FASTQ="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/blood_1C01_hifi_reads.fastq.gz"

# Input sample names
SAMPLE_NAME="tumor-all"
CONTROL_NAME="blood"

# Create necessary directories
mkdir -p "$PWD/fastq/"
mkdir -p "$PWD/bam/$SAMPLE_NAME/"
mkdir -p "$PWD/bam/$CONTROL_NAME/"
mkdir -p "$PWD/output/$SAMPLE_NAME/"
mkdir -p "$PWD/output/$CONTROL_NAME/"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PARSE_DIR"

# Step 1: Decompress FASTQ files
echo "Decompressing FASTQ files..."
gunzip -c "$SAMPLE_FASTQ" > "$PWD/fastq/$SAMPLE_NAME.fastq"
gunzip -c "$CONTROL_FASTQ" > "$PWD/fastq/$CONTROL_NAME.fastq"

# Step 2: Align sample reads to reference genome
echo "Aligning sample reads ($SAMPLE_NAME) with Minimap2 (PacBio HiFi mode)..."
minimap2 -ax map-hifi -t 32 "$REFERENCE" "$PWD/fastq/$SAMPLE_NAME.fastq" | \
    samtools view -Shb > "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam"

# Step 3: Sort and index BAM
echo "Sorting and indexing BAM for $SAMPLE_NAME..."
samtools sort -@ 32 -m 2G "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam" \
 -o "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.bam"
samtools index "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.bam"

# Step 4: Repeat for control sample
echo "Aligning control reads ($CONTROL_NAME) with Minimap2 (PacBio HiFi mode)..."
minimap2 -ax map-hifi -t 32 "$REFERENCE" "$PWD/fastq/$CONTROL_NAME.fastq" | \
    samtools view -Shb > "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam"

echo "Sorting and indexing BAM for $CONTROL_NAME..."
samtools sort -@ 32 -m 2G "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam" \
 -o "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.bam"
samtools index "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.bam"

# Step 5: Parse SVs using NanoMonsv
echo "Parsing structural variants for $SAMPLE_NAME..."
nanomonsv parse "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.bam" "$PWD/output/$SAMPLE_NAME/$SAMPLE_NAME"

echo "Parsing structural variants for $CONTROL_NAME..."
nanomonsv parse "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.bam" "$PWD/output/$CONTROL_NAME/$CONTROL_NAME"

# Step 6: Get SVs using NanoMonsv with control sample
echo "Running NanoMonsv SV detection..."
nanomonsv get \
    "$PWD/output/$SAMPLE_NAME/$SAMPLE_NAME" \
    "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.bam" \
    "$REFERENCE" \
    --control_prefix "$PWD/output/$CONTROL_NAME/$CONTROL_NAME" \
    --control_bam "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.bam" \
    --single_bnd \
    --use_racon \
    --processes 32 \
    --min_tumor_variant_read_num 2 \
    --min_tumor_VAF 0.01 \
    --max_control_variant_read_num 3 \
    --cluster_margin_size 100 \
    --median_mapQ_thres 25 \
    --max_overhang_size_thres 30 \
    --var_read_min_mapq 10

# Step 7: Cleanup
echo "Cleaning up temporary files..."
rm "$PWD/fastq/$SAMPLE_NAME.fastq"
rm "$PWD/fastq/$CONTROL_NAME.fastq"
rm "$PWD/bam/$SAMPLE_NAME/$SAMPLE_NAME.unsorted.bam"
rm "$PWD/bam/$CONTROL_NAME/$CONTROL_NAME.unsorted.bam"

echo "Nanomonsv SV calling completed."
