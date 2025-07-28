#!/bin/bash
### parameters required for teh LSF job ###
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J savana
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-savana.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-savana.err

### activate the savana conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate savana

### directories needed as well as reference and bam files ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/blood_1C01_hifi_reads_pbmm2.bam"
input_tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana"

### create the output directory ###
mkdir -p "$output_dir"

### call the variants with the homologous setup as in other callers ###
savana \
  --tumour "$input_tumor" \
  --normal "$input_blood" \
  --outdir "$output_dir" \
  --ref "$reference" \
  --pb \
  --mapq 10 \
  --min_support 2 \
  --min_af 0.01 \
  --length 20 \
  --buffer 100 \
  --insertion_buffer 100 \
  --single_bnd \
  --threads 32 \
  --predict_germline \
  --somatic_output "$output_dir/savana.somatic.vcf" \
  --germline_output "$output_dir/savana.germline.vcf"
