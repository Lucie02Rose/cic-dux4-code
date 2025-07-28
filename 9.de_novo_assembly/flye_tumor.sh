#!/bin/bash
### parameters for the LSF ###
#BSUB -n 64
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J flye_tumor_assembly
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-flye.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-flye.err

### activate the conda environment with flye ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate flye_env

### define the input, output ###
input_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_4_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/Flye"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/Flye/tmp"

### create the output and temporary directories ###
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### run flye with the pacbio preset, fastq, output directory, genome size and threads ###
flye --pacbio-raw "$input_fastq" --out-dir "$output_dir" --genome-size 3g --threads 64
### print statement to check outputs ###
echo "Flye assembly process has started. Check logs in $output_dir"
