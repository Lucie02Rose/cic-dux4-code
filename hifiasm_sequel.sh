#!/bin/bash
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sequeldenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sequeldenovo.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# List of input FASTQ.gz files
input_fastq1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A01_hifi_reads.fastq.gz"
input_fastq2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A02_hifi_reads.fastq.gz"
input_fastq3="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_2B01_hifi_reads.fastq.gz"
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_sequel_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-denovo"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-denovo/tmp"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Combine FASTQ.gz files
cat "$input_fastq1" "$input_fastq2" "$input_fastq3" > "$combined_fastq"

# Now run hifiasm with combined reads
hifiasm -o "$output_dir" -t64 "$combined_fastq"

echo "assembly process has started. Check logs in $output_dir"
