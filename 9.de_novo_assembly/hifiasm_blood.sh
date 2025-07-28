#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 64
#BSUB -M 150000
#BSUB -R 'span[hosts=1] select[mem>150000] rusage[mem=150000]'
#BSUB -q hugemem
#BSUB -J hifiblood
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-blooddenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-blooddenovo.err

### activate the hifiasm conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

### use the already 30x revio sample ###
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/blood_1C01_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/tmp"

### create the output and temporary directory, export the temporary directory ####
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### run hifiasm on the blood ###
hifiasm -o "$output_dir" -t64 "$combined_fastq"
### note that assembly is running ###
echo "assembly process has started"
