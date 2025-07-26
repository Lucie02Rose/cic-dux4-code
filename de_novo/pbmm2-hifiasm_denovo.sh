#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q long
#BSUB -J denovo
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
mom="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/mom_denovo/mom_denovo.fasta"
blood="/lustre/scratch126/casm/team274sb/lr26/a"
tumor=

output_mom="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/mom_denovo"
output_blood="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/blood_denovo"
tmp_dir_mom="$output_mom/tmp"
tmp_dir_blood="$output_blood/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir_mom"
mkdir -p "$tmp_dir_blood"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir_mom"  # Default set for mom

# Run pbmm2 Alignment for mom
echo "Aligning $mom to $reference..."
pbmm2 align "$reference" "$mom" "$output_mom/mom_aligned.bam" --preset CCS --sort -j 16 --unmapped

echo "pbmm2 alignment for mom completed: $output_mom/mom_aligned.bam"

# Now run for blood
export TMPDIR="$tmp_dir_blood"  # Change TMPDIR for blood

echo "Aligning $blood to $reference..."
pbmm2 align "$reference" "$blood" "$output_blood/blood_aligned.bam" --preset CCS --sort -j 16 --unmapped

echo "pbmm2 alignment for blood completed: $output_blood/blood_aligned.bam"

