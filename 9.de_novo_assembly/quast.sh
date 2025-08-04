#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J quast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-quast.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-quast.err

### activate the quast conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate quast_busco

### define directories fo all files to evaluate ###
flye="/lustre/scratch126/cellgen/behjati/lr26/Flye/assembly.fasta"
sequel="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-sequel.bp.p_ctg.fasta"
sequel_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-sequel.bp.hap1.p_ctg.fasta"
sequel_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/PacBio-sequel.bp.hap2.p_ctg.fasta"
revio="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/PacBio-revio.bp.p_ctg.fasta"
revio_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/PacBio-revio.bp.hap1.p_ctg.fasta"
revio_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/PacBio-revio.bp.hap2.p_ctg.fasta"
blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/PacBio-blood.bp.p_ctg.fasta"
blood_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/PacBio-blood.bp.hap1.p_ctg.fasta"
blood_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/PacBio-blood.bp.hap2.p_ctg.fasta"
mom="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/PacBio-mom.bp.p_ctg.fasta"
mom_1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/PacBio-mom.bp.hap1.p_ctg.fasta"
mom_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/PacBio-mom.bp.hap2.p_ctg.fasta"

### define and make the output directory ###
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-quast"

mkdir -p "$output_dir"

### run quast on all assemblies ###
quast -o "$output_dir" "$flye" "$sequel" "$sequel_1" "$sequel_2" "$revio" "$revio_1" "$revio_2" "$blood" "$blood_1" "$blood_2" "$mom" "$mom_1" "$mom_2"






