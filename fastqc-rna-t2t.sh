#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/fastqc%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/fastqc%J.e
#BSUB -n 8
#BSUB -M 20000
#BSUB -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]'
#BSUB -q normal
#BSUB -J fastqc-fastq
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate base

### Directories
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA/Fastqc"
input_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
rna_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"

mkdir -p "$output_dir"

find "$rna_dir" -name "*.fastq.gz" | \
    parallel -j 8 fastqc {} -o "$output_dir"



