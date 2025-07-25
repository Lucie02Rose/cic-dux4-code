#!/bin/bash
#BSUB -n 32
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-nanomonsv-match
#BSUB -G team274  # Fix extra space
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Activate conda environment for Nanomonsv
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate nanomonsv

# Define paths
PARSE_DIR="/lustre/scratch126/casm/team274sb/lr26/nanomonsv-parse"
OUTPUT_DIR="/lustre/scratch126/casm/team274sb/lr26/nanomonsv-results-matched"
REFERENCE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
SAMPLE_FASTQ="/lustre/scratch126/casm/team274sb/lr26/fastq-revio/tumor-all.fastq.gz"
CONTROL_FASTQ="/lustre/scratch126/casm/team274sb/lr26/fastq-revio/m84047_240202_152510_s2.hifi_reads.bc2025.fastq.gz"

# Input sample names
SAMPLE_NAME="tumor-all"
CONTROL_NAME="m84047_240202_152510_s2"

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


(... more ...)
------------------------------------------------------------

Successfully completed.

Resource usage summary:

    CPU time :                                   304976.00 sec.
    Max Memory :                                 65794 MB
    Average Memory :                             19372.94 MB
    Total Requested Memory :                     200000.00 MB
    Delta Memory :                               134206.00 MB
    Max Swap :                                   1 MB
    Max Processes :                              56
    Max Threads :                                93
    Run time :                                   28371 sec.
    Turnaround time :                            198478 sec.

The output (if any) is above this job summary.



PS:

Read file </lustre/scratch126/casm/team274sb/lr26/error_logs/316121.err> for stderr output of this job.

