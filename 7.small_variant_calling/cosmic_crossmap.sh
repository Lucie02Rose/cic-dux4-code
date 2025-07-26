#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J liftover
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/liftover.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/liftover.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
cosmic="/lustre/scratch126/casm/team274sb/lr26/T2T/cosmic_sorted_final_fixed_hg38.vcf"
chain="/lustre/scratch126/casm/team274sb/lr26/T2T/grch38-chm13v2.chain"
T2T="/lustre/scratch126/casm/team274sb/lr26/T2T"

cd "$T2T"

CrossMap vcf "$chain" \
    "$cosmic" \
    "$reference" \
    Cosmic_MutantCensus_T2T.vcf














