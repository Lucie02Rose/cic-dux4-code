#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J pbmm2-sniffles2-t2t-2B01-bam
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate sniffles2_env

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_bam_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment/2_B01"
output_vcf_dir="/lustre/scratch126/casm/team274sb/lr26/sniffles2-results/2_B01"

# Create output directory if it does not exist
mkdir -p "$output_vcf_dir"

# Process each BAM file in the directory
for bam_file in "$input_bam_dir"/*.bam; do
    # Extract base name without extension
    base_name=$(basename "$bam_file" .bam)
    output_vcf="$output_vcf_dir/${base_name}.vcf"

    echo "Calling SVs with Sniffles on $bam_file using $reference..."

    sniffles --input "$bam_file" \
             --vcf "$output_vcf" \
             --mosaic \
             --reference "$reference" \
             --threads 16 \
             --output-rnames

    echo "Sniffles SV calling completed: $output_vcf"
done
