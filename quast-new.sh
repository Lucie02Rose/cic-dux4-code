#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J quast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-quast.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-quast.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate quast_busco

### Directories
sequel="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-denovo.bp.p_ctg.fasta"
sequel_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-denovo.bp.hap1.p_ctg.fasta"
sequel_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-denovo.bp.hap2.p_ctg.fasta"
tumall_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor-whole/PacBio-tumor-whole.bp.hap1.p_ctg.fasta"
tumall_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor-whole/PacBio-tumor-whole.bp.hap2.p_ctg.fasta"
tumall="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor-whole/PacBio-tumor-whole.bp.p_ctg.fasta"

# Define the output directory
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-quast"

mkdir -p "$output_dir"

# Run QUAST on all assemblies
quast -o "$output_dir" "$sequel" "$sequel_1" "$sequel_2" "$tumall" "$tumall_1" "$tumall_2"






