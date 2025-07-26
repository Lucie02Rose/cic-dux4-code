#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/fastqc%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/fastqc%J.e
#BSUB -n 8
#BSUB -M 40000
#BSUB -R 'span[hosts=1] select[mem>40000] rusage[mem=40000]'
#BSUB -q normal
#BSUB -J fastqc-fastq
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate base

### Directories
output_dir="/lustre/scratch126/casm/team274sb/lr26/fastqc-rna"
input_dir="/lustre/scratch126/casm/team274sb/lr26/Wachtel/Wachtel/"
rna_dir="/lustre/scratch126/casm/team274sb/lr26/rna-t2t/"

mkdir -p "$output_dir"

find "$input_dir" -name "*.fastq.gz" | \
    parallel -j 8 fastqc {} -o "$output_dir"

# Ensure output directory exists
cd "$input_dir"


