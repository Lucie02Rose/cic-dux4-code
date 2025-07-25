#!/bin/bash
#BSUB -n 32
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-aligtumor.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-aligtumor.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# List of input FASTQ.gz files
reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"
#input_fastq1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_sequel_hifi_reads.fastq.gz"
#input_fastq2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"

combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_4_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/tmp"

# Combine FASTQ.gz files
#cat "$input_fastq1" "$input_fastq2" > "$combined_fastq"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir"

# Define output BAM file name
base_name=$(basename "$combined_fastq" .fastq.gz)
output_bam="$output_dir/${base_name}_pbmm2.bam"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# Run pbmm2 Alignment
echo "Aligning $input_fastq to $reference..."
pbmm2 align "$reference" "$combined_fastq" "$output_bam" --preset HIFI --sort -j 32 --unmapped

echo "pbmm2 alignment completed: $output_bam"


echo "assembly process has started. Check logs in $output_dir"
