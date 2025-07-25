#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-alignment-revio-fastq
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
input_dir="/lustre/scratch126/casm/team274sb/lr26/fastq-revio"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/revio"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$output_dir" "$tmp_dir" 

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# run iteratively for all files
for input_fastq in "$input_dir"/*fastq.gz; do

        base_name=$(basename "$input_fastq".fastq.gz)
        output_bam="$output_dir/${base_name}_pbmm2-farm22.bam"

        # Run the pbmm2 alignment

        echo "Aligning $input_fastq to $reference..."

        pbmm2 align "$reference" "$input_fastq" "$output_bam" --preset CCS --sort ->

        echo "pbmm2 alignment completed: $output_bam"

        done

