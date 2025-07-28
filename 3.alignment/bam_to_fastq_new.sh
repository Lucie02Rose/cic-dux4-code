#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J bam2fastq
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/bamfastq-%J.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/bamfastq-%J.e

### activate the conda environment with the basm2fastq in it ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### input and output directory ###
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq"
### defining sample names ###
samples=("blood_1C01_hifi_reads.bam" "tumor_1A01_hifi_reads.bam" "tumor_1A02_hifi_reads.bam" "tumor_1B01_hifi_reads.bam" "tumor_2B01_hifi_reads.bam" "mom_1B02_hifi_reads.bam")
### make the output directory ###
mkdir -p "$output_dir"
### for each of the samples define the bam file path and name ###
for bam_file_name in "${samples[@]}"; do
  bam_file="$input_dir/$bam_file_name"
  ### check if they are there ###
  if [ -f "$bam_file" ]; then
    echo "Processing file: $bam_file"
    base_name=$(basename "$bam_file" .bam)

    ### convert bam to fastq ###
    bam2fastq "$bam_file" -o "$output_dir/$base_name"

    ### check if the output exists and then compress ###
    if [ -f "$output_dir/$base_name.fastq" ]; then
      if command -v pigz &> /dev/null; then
        pigz "$output_dir/$base_name.fastq"
      else
        gzip "$output_dir/$base_name.fastq"
      fi
    else
    ### print statement ###
      echo "Warning: FASTQ file for $bam_file was not created."
    fi
  else
  ### print statement ###
    echo "File $bam_file does not exist. Skipping."
  fi
done
