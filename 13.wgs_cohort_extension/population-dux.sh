#!/bin/bash
#BSUB -n 40
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J t2t-alignment-single
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/single.%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/single.%J.e

set -euo pipefail

echo "Activating environment"
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bwa-mem2

# Set parameters
threads=40
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
roi_region="chr4:193541548-193563474"
roi_fasta="chr4_dux_augustus.fasta"
output_dir="/lustre/scratch126/casm/team274sb/lr26/population_dux"
input_bam="/lustre/scratch126/casm/team274sb/lr26/population_dux/PD54859b/PD54859b.sample.dupmarked.bam"

# Prepare workspace
sample=$(basename "$input_bam" .bam)
cd "$output_dir"   # Extract ROI fasta once
if [[ ! -f "$roi_fasta" ]]; then
    echo "Extracting ROI from reference"
    samtools faidx -@ "$threads" "$reference" "$roi_region" > "$roi_fasta"
fi

echo "Processing sample: $sample"

samtools fastq -@ "$threads" -1 "${sample}_R1.fq" -2 "${sample}_R2.fq" "$input_bam"

# Align reads to ROI
minimap2 -t "$threads" -x sr "$roi_fasta" "${sample}_R1.fq" "${sample}_R2.fq" > "${sample}_roi.paf"

# Filter for matching reads
awk '$10/$2 >= 0.25 {print $1}' "${sample}_roi.paf" | sort | uniq > "${sample}_passing_ids.txt"
if [[ ! -s "${sample}_passing_ids.txt" ]]; then
    echo "No reads passed for $sample — skipping."
    exit 0
fi

# Subsample BAM
samtools view -@ "$threads" -N "${sample}_passing_ids.txt" -b "$input_bam" > "${sample}_temp.bam"
samtools sort -@ "$threads" -n "${sample}_temp.bam" -o "${sample}_temp.sorted.bam"
samtools fixmate -@ "$threads" -m "${sample}_temp.sorted.bam" "${sample}_fixmate.bam"
samtools view -@ "$threads" -bf 0x2 "${sample}_fixmate.bam" > "${sample}_subsampled.bam"  # FASTQ from subsample
samtools fastq -@ "$threads" -1 "${sample}_sub_R1.fq" -2 "${sample}_sub_R2.fq" "${sample}_subsampled.bam"

# Final alignment to full reference
bwa-mem2 mem -t "$threads" "$reference" "${sample}_sub_R1.fq" "${sample}_sub_R2.fq" | \
    samtools view -@ "$threads" -Sb - | \
    samtools sort -@ "$threads" -o "${sample}_t2t.sorted.bam"
samtools index -@ "$threads" "${sample}_t2t.sorted.bam"

# Cleanup
rm -f "${sample}_R1.fq" "${sample}_R2.fq" "${sample}_roi.paf" "${sample}_passing_ids.txt"
rm -f "${sample}_temp.bam" "${sample}_temp.sorted.bam" "${sample}_fixmate.bam" "${sample}_subsampled.bam"
rm -f "${sample}_sub_R1.fq" "${sample}_sub_R2.fq"
