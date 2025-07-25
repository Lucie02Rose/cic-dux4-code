#!/bin/bash
#BSUB -n 16
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q basement
#BSUB -J pbmm2-alignment-t2t-tumor-all
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
echo "Activating bioinfo environment"
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

# Directories
echo "Handling directories"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
s1A01="/lustre/scratch126/casm/team274sb/lr26/PacBio_raw/1_A01/m64094e_230126_154129.hifi_reads.bam"
s1A02="/lustre/scratch126/casm/team274sb/lr26/PacBio_raw/1_A02/m64178e_230206_134948.hifi_reads.bam"
s2B01="/lustre/scratch126/casm/team274sb/lr26/PacBio_raw/2_B01/m64178e_230207_165902.hifi_reads.bam"
s1B01="/lustre/scratch126/casm/team274sb/lr26/Revio_raw/1_B01/m84047_230404_172053_s2.hifi_reads.default.bam"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-tumor-all-from-bams"
tmp_dir="$output_dir/tmp"

# Create temp directory if it doesn't exist
echo "Creating directories"
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"

# Define output BAM file name
echo "Merging bam files"
samtools merge -@ 16 -u "$output_dir/tumor_merged_raw.bam" "$s1A01" "$s1A02" "$s2B01" "$s1B01"

echo "Creating output file name"
output_bam="$output_dir/tumor_final_merged_wtags.bam"

# Set TMPDIR for pbmm2 (samtools sorting)
echo "Setting tmp dir"
export TMPDIR="$tmp_dir"

echo "Activating new env with pbmm2"
conda deactivate
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Run pbmm2 Alignment
echo "Aligning to $reference..."
pbmm2 align "$reference" "$output_dir/tumor_merged_raw.bam" "$output_bam" --preset CCS --sort -j 16 --unmapped

echo "pbmm2 alignment completed: $output_bam"

