#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span [hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal 
#BUSB -J indexing

### activate my conda environment - pbcore is a separate and used for detailed metrics and qc
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference_gz="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa.gz"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/lr26/T2T"
output_dir="/lustre/scratch126/casm/team274sb/lr26/T2T"

# Decompress reference if not already done
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi
