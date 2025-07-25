#!/bin/bash
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/hg004.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/hg004.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# Define directories
input_fastq="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/fastq-revio-sprq/m84039_241002_020632_s4.hifi_reads.bc2021.fastq.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_denovo_hifiasm"
tmp_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_denovo_hifiasm/denovhg004"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Run Flye with reduced memory & threads
hifiasm -o "$output_dir" -t64 "$input_fastq"

echo "Flye assembly process has started. Check logs in $output_dir"

