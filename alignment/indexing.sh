#!/bin/bash
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-index.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-index.e
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BUSB -J indexing
#BSUB -G team274

### activate my conda environment - pbcore is a separate and used for detailed metrics and qc
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference_gz="/lustre/scratch125/cellgen/behjati/lr26/T2T/chm13v2.0.fa.gz"
reference="/lustre/scratch125/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
ref_dir="/lustre/scratch125/cellgen/behjati/lr26/T2T"

cd "$ref_dir"

# Decompress reference if not already done
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi

# Indexing the T2T reference genome (if index files don't exist)
if [ ! -f "${reference%.fa}.mmi" ]; then
    pbmm2 index "$reference" "${reference%.fa}.mmi"
fi
