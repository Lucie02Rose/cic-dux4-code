#!/bin/bash
#BSUB -n 16
#BSUB -M 10000
#BSUB -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]'
#BSUB -q normal
#BSUB -J compression
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
gzip /lustre/scratch126/casm/team274sb/lr26/pbmm2-tumor-all-from-bams/tumor_final_merged_wtags.bam

echo "gzip completed"
