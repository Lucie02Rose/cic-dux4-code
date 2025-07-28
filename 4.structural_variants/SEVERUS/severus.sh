#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J severus
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-severus.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-severus.err

### activate the severus conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate severus_env

### define the input, outpu and reference ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
input_blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/blood_1C01_hifi_reads_pbmm2.bam"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus"
### create the output directory ###
mkdir -p "$output_dir"
### print statement for processing ###
echo "Running Severus..."
### same settings as with the previous somatic callers ###
severus \
   --target-bam "$input_tumor" \
   --control-bam "$input_blood" \
   --out-dir "$output_dir" \
   -t 32 \
   --min-support 2 \
   --vaf-thr 0.01 \
   --min-mapq 10 \
   --bp-cluster-size 100 \
   --min-sv-size 20 \
   --single-bp \
   --output-read-ids
### pring statement for being finished ###
echo "Severus SV call completed"
