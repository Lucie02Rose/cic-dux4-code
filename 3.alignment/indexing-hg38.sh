#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J indexing-hg38
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-index-hg38.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-index-hg38.e

### activate my conda environment base with pbmm2 ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### directories ###
reference_gz="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.fa.gz"
reference="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.fa"
ref_dir="/lustre/scratch126/cellgen/behjati/lr26/hg38"
### change to reference ###
cd "$ref_dir"

### possibly decompress reference ###
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi

### index the reference if not already done ###
if [ ! -f "${reference%.fa}.mmi" ]; then
    pbmm2 index "$reference" "${reference%.fa}.mmi"
fi

