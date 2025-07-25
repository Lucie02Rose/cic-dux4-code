#!/bin/bash
#BSUB -n 32
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J pepper
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Load Singularity module
module load ISG/singularity/3.11.4

# Set cache and temp directories (modify paths as needed)
export SINGULARITY_CACHEDIR=/lustre/scratch126/casm/team274sb/lr26/singularity
export SINGULARITY_TMPDIR=/lustre/scratch126/casm/team274sb/lr26/singularity/tmp

# Define paths
SINGULARITY_IMG="/lustre/scratch126/casm/team274sb/lr26/singularity/pepper_deepvariant_r0.8.sif"
BAM_FILE="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/1_C01-revio/m84047_240202_155616_s3.hifi_reads.bc2026_pbmm2.bam"
REFERENCE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
OUTPUT_DIR="/lustre/scratch126/casm/team274sb/lr26/pepper-blood"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Run PEPPER-Margin-DeepVariant with Singularity
singularity exec --bind /lustre/scratch126/casm/team274sb/lr26:/mnt \
    "${SINGULARITY_IMG}" \
    run_pepper_margin_deepvariant call_variant \
    -b "/mnt/pbmm2-alignment/1_C01-revio/m84047_240202_155616_s3.hifi_reads.bc2026_pbmm2.bam" \
    -f "/mnt/T2T/chm13v2.0.fa" \
    -o "/mnt/pepper-blood" \
    -p "blood_output" \
    -t 32 \
    --hifi
