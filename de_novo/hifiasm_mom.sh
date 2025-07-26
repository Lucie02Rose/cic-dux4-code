#!/bin/bash
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-momdenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-momdenovo.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# List of input FASTQ.gz files
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/mom_1B02_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/tmp"

# Create temporary directory
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Now run hifiasm with combined reads
hifiasm -o "$output_dir" -t64 "$combined_fastq"

echo "assembly process has started. Check logs in $output_dir"
