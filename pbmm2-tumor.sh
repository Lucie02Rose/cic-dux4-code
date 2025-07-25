#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q long
#BSUB -J pbmm2-combined-tumor-fastq-t2t
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J-%J_NAME.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J-%J_NAME.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
input_dir="/lustre/scratch126/casm/team274sb/lr26/fastq-revio/tumor-all.fastq.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-tumor"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$output_dir" "$tmp_dir" 

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir"

# run iteratively for all files
for input_fastq in "$input_dir"/*fastq.gz; do

        base_name=$(basename "$input_fastq".fastq.gz)
        output_bam="$output_dir/${base_name}_pbmm2-tumor-allt2t.bam"

        # Run the pbmm2 alignment

        echo "Aligning $input_fastq to $reference..."

        pbmm2 align "$reference" "$input_fastq" "$output_bam" --preset CCS --sort -j 16 --unmapped

        echo "pbmm2 alignment completed: $output_bam"

        done

