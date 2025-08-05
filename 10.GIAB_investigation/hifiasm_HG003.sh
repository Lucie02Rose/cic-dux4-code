#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 64
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q hugemem
#BSUB -J hg003
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-hg003.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-hg003.err

### activate the hifiasm conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# Define directories
input_fastq="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003-rep1.fastq"
output_dir="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_revio"
tmp_dir="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_revio/denovhg003"

### create the output and temporary directory, export the temporary directory ####
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### run hifiasm ###
hifiasm -o "$output_dir" -t64 "$input_fastq"
### note that assembly is running ###
echo "Assembly process has started. Check logs in $output_dir"
