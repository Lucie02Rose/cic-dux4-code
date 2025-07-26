#!/bin/bash
#BSUB -n 16
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J bloodhif
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
blood="/lustre/scratch126/casm/team274sb/lr26/hifiasm_blood_denovo"
output_blood="/lustre/scratch126/casm/team274sb/lr26/hifiasm_blood_denovo"
tmp_dir_blood="$output_blood/tmp"

# Create temp directory if it doesn't exist
mkdir -p "$tmp_dir_blood"

# Set TMPDIR for pbmm2 (samtools sorting)
export TMPDIR="$tmp_dir_blood"  # Default set for mom

# Run pbmm2 Alignment for mom
echo "Aligning  to $reference..."

# Loop through all BAM files in the input directory and subdirectories
for input_fasta in "$blood"/*p_ctg.fasta; do
    # Extract sample name from BAM filename
    sample_name=$(basename "$input_fasta" .fasta)

    echo "Processing: $input_fasta"

    # Run methylation analysis
    pbmm2 align "$reference" "$input_fasta" "$output_blood/${sample_name}_aligned.bam" --preset CCS --sort -j 16 --unmapped 

    echo "Completed: $input_fasta"
done

echo "pbmm2 alignments for blood completed."


