#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span [hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal 
#BSUB -J pbmm2indexing-hg38
#BSUB -G team274

### activate my conda environment - pbcore is a separate and used for detailed metrics and qc
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
input_dir="/lustre/scratch126/casm/team274sb/lr26/hg38"
output_dir="/lustre/scratch126/casm/team274sb/lr26/hg38"

# Indexing the T2T reference genome (if index files don't exist)
pbmm2 index "$reference" "${reference%.fa}.mmi"
