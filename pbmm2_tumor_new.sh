#!/bin/bash
#BSUB -n 32
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q long
#BSUB -J tumor-all
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-tumor.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-tumor.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"
input_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir"

# Define output BAM file name
base_name=$(basename "$input_fastq" .fastq.gz)
output_bam="$output_dir/${base_name}_pbmm2.bam"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Run pbmm2 Alignment
echo "Aligning $input_fastq to $reference..."
pbmm2 align "$reference" "$input_fastq" "$output_bam" --preset HIFI --sort -j 32 --unmapped

echo "pbmm2 alignment completed: $output_bam"
