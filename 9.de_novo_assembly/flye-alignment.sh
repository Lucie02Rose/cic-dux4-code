#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J flye-alignment
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-flye-alignment.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-flye-alignment.e

### activate the environment with minimap and samtools ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base
### references, input and output ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.mmi"
input="/lustre/scratch126/cellgen/behjati/lr26/Flye/assembly.fasta"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/Flye"

### alignment with minimap ###
echo "Aligning assembly to t2t..."
minimap2 -ax asm20 "$reference" "$input" -o "$output_dir/all_contigs_vs_T2T.sam"

### alignment completed and then sort and index the sam file ###
echo "Alignment completed, sorting and indexing with samtools"
### convert to bam file ###
samtools view -Sb "$output_dir/all_contigs_vs_T2T.sam" > "$output_dir/all_contigs_vs_T2T.bam"
### sort ###
samtools sort "$output_dir/all_contigs_vs_T2T.bam" -o "$output_dir/all_contigs_vs_T2T_sorted.bam"
### index ###
samtools index "$output_dir/all_contigs_vs_T2T_sorted.bam"
### message that it is completed ###
echo "Completed."
