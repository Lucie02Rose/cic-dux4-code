#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J sawfish-hg38-tumor
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J-sawfish-tumor-hg38.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J-sawfish-tumor-hg38.err

### activate the sawfish conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate sawfish

#### define reference, input and output directories ###
bam_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned-hg38/tumor_all_4_hifi_reads_pbmm2.bam"
reference="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.fa"
output_vcf_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-hg38"

### extract the basename from the file ###
base_name=$(basename "$bam_file" .bam)
### define the two output directories needed by sawfish ###
discover_dir="$output_vcf_dir/${base_name}_discover"
joint_call_dir="$output_vcf_dir/${base_name}_joint_call"
### remove any previous ones for this file ###
[ -d "$discover_dir" ] && rm -rf "$discover_dir"
[ -d "$joint_call_dir" ] && rm -rf "$joint_call_dir"
### print statement to run the discover step ###
echo "Running Sawfish discover on $bam_file using $reference..."
### discover structural variants step ###
sawfish discover --bam "$bam_file" --output-dir "$discover_dir" --ref "$reference" --threads 16
### print statement to run the second joint call step ###
echo "Running Sawfish joint-call for $base_name..."
### run the second joint call step ###
sawfish joint-call --sample "$discover_dir" --output-dir "$joint_call_dir" --threads 16
### print statement for finished ###
echo "Sawfish processing completed for: $base_name"
