#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q normal
#BSUB -J liftover
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/liftoverv.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/liftoverv.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
cosmic="/lustre/scratch126/casm/team274sb/lr26/T2T/cosmic_sorted_final_fixed_hg38.vcf"
chain="/lustre/scratch126/casm/team274sb/lr26/T2T/grch38-chm13v2.chain"
T2T="/lustre/scratch126/casm/team274sb/lr26/T2T"

cd "$T2T"

picard -Xmx100g CreateSequenceDictionary \
    R="$reference" \
    O=chm13v2.0.dict

picard -Xmx100g LiftoverVcf \
    I="$cosmic" \
    O=cosmic_t2t.vcf \
    CHAIN="$chain" \
    REJECT=unlifted_cosmic.vcf \
    R="$reference"













