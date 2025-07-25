#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q long
#BSUB -J pbmm2-alignment-t2t-lucie-1A01-new-bam
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
input_fastq="/lustre/scratch126/casm/team274sb/lr26/PacBio_raw/1_A01/m64094e_230126_154129.hifi_reads.bam"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/1_A01"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir"

# Define output BAM file name
base_name=$(basename "$input_fastq" .hifi_reads.fastq.fastq.gz)
output_bam="$output_dir/${base_name}_pbmm2-farm22-bam.bam"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Run pbmm2 Alignment
echo "Aligning $input_fastq to $reference..."
pbmm2 align "$reference" "$input_fastq" "$output_bam" --preset CCS --sort -j 16 --unmapped

echo "pbmm2 alignment completed: $output_bam"
