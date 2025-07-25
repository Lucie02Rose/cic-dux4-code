#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J pepper-mom
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/pepper-intersec.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/pepper-intersec.err

source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

# Define pathways
mom_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-mom/mom_output.vcf.gz"
blood_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-blood/blood_output.vcf.gz"
tumor_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/tumor_output.vcf.gz"
mom_norm="/lustre/scratch126/casm/team274sb/lr26/pepper-mom/mom_output_norm.vcf.gz"
blood_norm="/lustre/scratch126/casm/team274sb/lr26/pepper-blood/blood_output_norm.vcf.gz"
tumor_norm="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/tumor_output_norm.vcf.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01"

# Change to pepper
cd "$output_dir"

bcftools norm -m -any -Oz -o "$mom_norm" "$mom_vcf"
bcftools norm -m -any -Oz -o "$blood_norm" "$blood_vcf"
bcftools norm -m -any -Oz -o "$tumor_norm" "$tumor_vcf"
bcftools index "$mom_norm"
bcftools index "$blood_norm"
bcftools index "$tumor_norm"

# Unique to file1 (not in file2 or file3)
bcftools isec -p isec1 -C "$tumor_norm" "$blood_norm" "$mom_norm"

# Common to files2 and 3
bcftools isec -p isec_2_3_common -n=2 "$blood_norm" "$mom_norm"
