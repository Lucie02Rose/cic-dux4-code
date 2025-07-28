#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J 1B02
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-1B02.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-1B02.err

### activate the base conda environment with pbmm2 ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define the directories to be used ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_bam="/lustre/scratch126/cellgen/behjati/lr26/PacBio/mom_1B02_hifi_reads.bam"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="$output_dir/tmp"

### create directories which may not exist ###
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"

### output bam name ###
base_name=$(basename "$input_bam" .bam)
output_bam="$output_dir/${base_name}_pbmm2.bam"

### export the tmp dir for pbmm2 sorting ###
export TMPDIR="$tmp_dir"

### run the pbmm2 alignment with sorting, HIFI and multithreaded and unmapped flags ###
echo "Aligning $input_bam to $reference..."
pbmm2 align "$reference" "$input_bam" "$output_bam" --preset HIFI --sort -j 16 --unmapped
### end process message ###
echo "pbmm2 alignment completed: $output_bam"
