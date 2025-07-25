#!/bin/bash

#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 4
#BSUB -M 4000
#BSUB -R 'span[hosts=1] select[mem>4000] rusage[mem=4000]'
#BSUB -q normal
#BSUB -J fastqc-fastq-all
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate base

### Directories
input_dir="/lustre/scratch126/casm/team274sb/lr26/fastq"

# Ensure output directory exists
cd "$input_dir"

zcat m64094e_230126_154129.hifi_reads.fastq.fastq.gz m64178e_230206_134948.hifi_reads.fastq.fastq.gz m64178e_230207_165902.hifi_reads.fastq.fastq.gz | gzip > sequel-all.fastq.gz
