#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J indexing-t2t
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-index-t2t.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-index-t2t.e

### activate my conda environment - base has pbmm2 ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define the directories with the reference
reference_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa.gz"
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
ref_dir="/lustre/scratch126/cellgen/behjati/lr26/T2T"
### change to the reference directory ###
cd "$ref_dir"

### decompress the reference if it is not yet done ###
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi

### index the genome to generate .mmi if not already done ###
if [ ! -f "${reference%.fa}.mmi" ]; then
    pbmm2 index "$reference" "${reference%.fa}.mmi"
fi
