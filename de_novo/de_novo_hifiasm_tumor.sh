#!/bin/bash
#BSUB -n 64
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-denovotum.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-denovotum.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# Define directories
input_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-1B01denovo"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/denovtum"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Run Flye with reduced memory & threads
hifiasm -o "$output_dir" -t64 "$input_fastq"

echo "Flye assembly process has started. Check logs in $output_dir"

