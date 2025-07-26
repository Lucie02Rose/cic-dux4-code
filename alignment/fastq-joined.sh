#!/bin/bash

#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 4
#BSUB -M 40
#BSUB -R 'span[hosts=1] select[mem>40] rusage[mem=40]'
#BSUB -q normal
#BSUB -J fastqc-fastq-all
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate base

### Directories
input_dir="/lustre/scratch126/casm/team274sb/lr26/fastq-revio"

# Ensure output directory exists
cd "$input_dir"

zcat sequel-all.fastq.gz 1B01-all.fastq.gz | gzip > tumor-all.fastq.gz
