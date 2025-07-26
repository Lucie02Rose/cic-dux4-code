#!/bin/bash
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/hg002.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/hg002.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# Define directories
input_fasta="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/combined_HG003.fasta"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_new"
tmp_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_new/denovhg003"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
cd "$output_dir"
export TMPDIR="$tmp_dir"

# Run Flye with reduced memory & threads
hifiasm -o "$output_dir" -t64 "$input_fasta"

echo "Flye assembly process has started. Check logs in $output_dir"

