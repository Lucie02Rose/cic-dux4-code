#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 64
#BSUB -M 150000
#BSUB -R 'span[hosts=1] select[mem>150000] rusage[mem=150000]'
#BSUB -q hugemem
#BSUB -J tumor_hifi_revio
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-tumor_hifi_revio.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-tumor_hifi_revio.err

### activate the hifiasm conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

### use the already 30x revio sample ###
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/tmp"

### create the output and temporary directory, export the temporary directory ###
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### run hifiasm on the revio ###
hifiasm -o "$output_dir" -t64 "$combined_fastq"
### note that assembly is running ###
echo "assembly process has started"
