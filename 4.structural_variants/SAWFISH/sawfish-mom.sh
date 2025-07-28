#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J sawfish-mom
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sawfish-mom.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sawfish-mom.err

### activate the sawfish conda environment ###
source /software/cellgen/team274/lr26/miniforge3/etc/profile.d/conda.sh
conda activate sawfish
### define the reference, input and output ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
bam_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/mom_2B01_hifi_reads_pbmm2.bam"
output_vcf_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish"

### extract the sample name from the file ###
base_name=$(basename "$bam_file" .bam)
### define outputs for the directories that need to be created for sawfish ###
discover_dir="$output_vcf_dir/${base_name}_discover"
joint_call_dir="$output_vcf_dir/${base_name}_joint_call"
### remove any previous for this sample ###
[ -d "$discover_dir" ] && rm -rf "$discover_dir"
[ -d "$joint_call_dir" ] && rm -rf "$joint_call_dir"
### print statement that it is running ###
echo "Running Sawfish discover on $bam_file using $reference..."
### discover variants ###
sawfish discover --bam "$bam_file" --output-dir "$discover_dir" --ref "$reference" --threads 16
### print statement for the second step ###
echo "Running Sawfish joint-call for $base_name..."
### joint call from the discovery ###
sawfish joint-call --sample "$discover_dir" --output-dir "$joint_call_dir" --threads 16
### print statement for finished ###
echo "Sawfish processing completed for: $base_name"

