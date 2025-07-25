#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-assembly-t2t
#BSUB -G team274

source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
input="/lustre/scratch126/casm/team274sb/lr26/flye-pacbio-tumor/assembly.fasta"
output_dir="/lustre/scratch126/casm/team274sb/lr26/flye-pacbio-tumor/"

echo "Aligning assembly to t2t..."

minimap2 -ax asm20 "$reference" "$input" -o "$output_dir/all_contigs_vs_T2T.sam"

echo "Alignment completed, activating bioinfo"

source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

echo "Sorting and indexing with samtools"

samtools view -Sb "$output_dir/all_contigs_vs_T2T.sam" > "$output_dir/all_contigs_vs_T2T.bam"

samtools sort "$output_dir/all_contigs_vs_T2T.bam" -o "$output_dir/all_contigs_vs_T2T_sorted.bam"

samtools index "$output_dir/all_contigs_vs_T2T_sorted.bam"

echo "Completed."
