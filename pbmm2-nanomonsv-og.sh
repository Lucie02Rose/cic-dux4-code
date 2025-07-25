#!/bin/bash
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q normal
#BSUB -J pbmm2-sniffles2-t2t-1C01-both
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Activate conda environment for Nanomonsv
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate nanomonsv

echo "Running Nanomonsv..."

# Define paths
PARSE_DIR="/lustre/scratch126/casm/team274sb/lr26/nanomonsv-parse"
OUTPUT_DIR="/lustre/scratch126/casm/team274sb/lr26/nanomonsv-results-matched"
TUMOR_BAM="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/1_B01-revio/m84047_230404_172053_s2.hifi_reads.default_pbmm2-defaultbam.bam"
CONTROL_BAM="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/1_B02-revio/m84047_240202_152510_s2.hifi_reads.bc2025_pbmm2-defaultbam.bam"
REFERENCE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"

# Define output prefixes for parsed data
TUMOR_PREFIX="$PARSE_DIR/tumor"
CONTROL_PREFIX="$PARSE_DIR/control"

# Create required directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PARSE_DIR"

# Step 1: Parse Tumor BAM
echo "Parsing tumor BAM..."
nanomonsv parse "$TUMOR_BAM" "$TUMOR_PREFIX" --debug > parse_log.txt 2>&1

# Step 2: Parse Control BAM (if using matched normal)
echo "Parsing control BAM..."
nanomonsv parse "$CONTROL_BAM" "$CONTROL_PREFIX" --debug > parse_log.txt 2>&1

# Step 3: Run Nanomonsv SV Calling
echo "Running Nanomonsv get..."
nanomonsv get \
    --min_tumor_variant_read_num 2 \
    --min_tumor_VAF 0.01 \
    --max_control_variant_read_num 3 \
    --cluster_margin_size 100 \
    --median_mapQ_thres 25 \
    --max_overhang_size_thres 30 \
    --var_read_min_mapq 10 \
    --single_bnd \
    --use_racon \
    "$TUMOR_PREFIX" "$TUMOR_BAM" "$REFERENCE" \
    --control_prefix "$CONTROL_PREFIX" --control_bam "$CONTROL_BAM"

if [[ $? -eq 0 ]]; then
    echo "✅ Nanomonsv SV calling completed successfully!"
else
    echo "❌ ERROR: Nanomonsv get failed."
    exit 1
fi
echo "Nanomonsv SV calling completed. Results are in: $OUTPUT_DIR"
