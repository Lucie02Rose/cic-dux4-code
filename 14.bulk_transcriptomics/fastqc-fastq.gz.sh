#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 8
#BSUB -M 20000
#BSUB -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]'
#BSUB -q normal
#BSUB -J fastqc-fastq
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-fastqc.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-fastqc.e

### activate the base conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define input and output ###
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA/Fastqc"
rna_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"

### make the output directory ###
mkdir -p "$output_dir" 

### run fastqc in parallel ###
find "$rna_dir" -name "*.fastq.gz" | \
    parallel -j 8 fastqc {} -o "$output_dir"
