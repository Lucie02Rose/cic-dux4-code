#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J bam2fastq_job
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/giabbamfastq.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/giabbamfastq.e

### activate conda environment - bam2fastq is in base ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define input and output directories - input is external downloads ###
input_dir="/lustre/scratch126/casm/team274sb/ExternalData/byUser/lr26/PacBio_GIAB_trio/bams/"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/fastq-revio-sprq"

### define the specific sample names to use ###
### these are the new chemistry samples - sprq ###
samples=("m84039_241001_220042_s2.hifi_reads.bc2018.bam" "m84039_241002_000337_s3.hifi_reads.bc2020.bam" "m84039_241002_020632_s4.hifi_reads.bc2021.bam")
### make the output directory ###
mkdir -p "$output_dir"
### process each in a for loop saving its path, name ###
for bam_file_name in "${samples[@]}"; do
  bam_file="$input_dir/$bam_file_name"

  if [ -f "$bam_file" ]; then
    echo "Processing file: $bam_file"
    base_name=$(basename "$bam_file" .bam)

    ### convert bam to fastq ###
    bam2fastq "$bam_file" -o "$output_dir/$base_name"

    ### compress fastq ###
    if [ -f "$output_dir/$base_name.fastq" ]; then
      if command -v pigz &> /dev/null; then
        pigz "$output_dir/$base_name.fastq"
      else
        gzip "$output_dir/$base_name.fastq"
      fi
    else
    ### check if it was created or no ###
      echo "Warning: FASTQ file for $bam_file was not created."
    fi
  else
  ### check if the bam exists or no ###
    echo "File $bam_file does not exist. Skipping."
  fi
done
