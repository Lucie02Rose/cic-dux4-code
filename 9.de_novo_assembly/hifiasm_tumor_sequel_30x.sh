#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 64
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q long
#BSUB -J sequel_hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sequeldenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sequeldenovo.err

### Acactivate the hifiasm conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

### concatenate and list the fastq files (this may already have been done in the alignment step) ###
input_fastq1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A01_hifi_reads.fastq.gz"
input_fastq2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A02_hifi_reads.fastq.gz"
input_fastq3="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_2B01_hifi_reads.fastq.gz"
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_sequel_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sequel/tmp"

### create the output and temporary directory, export the temporary directory ####
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### concatenate the sequel II runs to have one 30x depth of coverage file ###
cat "$input_fastq1" "$input_fastq2" "$input_fastq3" > "$combined_fastq"
### run hifiasm on the sequel II ###
hifiasm -o "$output_dir" -t64 "$combined_fastq"
### note that assembly is running ###
echo "assembly process has started"
