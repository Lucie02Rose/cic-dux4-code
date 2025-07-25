#!/bin/bash
#BSUB -n 64
#BSUB -M 600000
#BSUB -R 'span[hosts=1] select[mem>600000] rusage[mem=600000]'
#BSUB -q hugemem
#BSUB -J flye_mom_assembly
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate flye_env

# Define directories
input_fastq="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/mom_unique.fasta"
output_dir="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/mom_denovo_chr1410"
tmp_dir="/lustre/scratch126/casm/team274sb/lr26/denovmom"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Run Flye with reduced memory & threads
flye --pacbio-raw "$input_fastq" --out-dir "$output_dir" --genome-size 3g --threads 64

echo "Flye assembly process has started. Check logs in $output_dir"

