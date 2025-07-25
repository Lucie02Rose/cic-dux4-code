#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J cosmic
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/cos38%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/cos38%J.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
cosmic="/lustre/scratch126/casm/team274sb/lr26/hg38/cleaned_cosmic_hg38.vcf.gz"
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
mom="/lustre/scratch126/casm/team274sb/lr26/pepper-mom38/mom_output.vcf.gz"
blood="/lustre/scratch126/casm/team274sb/lr26/pepper-blood38/blood_output.vcf.gz"
tumor="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor38/tumor_output.vcf.gz"
dir="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor38"

cd "$dir"

bcftools norm -f "$reference" -Oz -o normalized_tumor_output.vcf.gz "$tumor"
bcftools norm -f "$reference" -Oz -o normalized_blood_output.vcf.gz "$blood"
bcftools norm -f "$reference" -Oz -o normalized_mom_output.vcf.gz "$mom"
bcftools norm -f "$reference" -Oz -o normalized_cosmic_hg38.vcf.gz "$cosmic"

tabix -p vcf normalized_tumor_output.vcf.gz
tabix -p vcf normalized_blood_output.vcf.gz
tabix -p vcf normalized_mom_output.vcf.gz
tabix -p vcf normalized_cosmic_hg38.vcf.gz

bcftools isec -p isec_tumor_germline -C normalized_tumor_output.vcf.gz normalized_blood_output.vcf.gz normalized_mom_output.vcf.gz
bcftools isec -p isec_blood_mom_somatic -n=2 normalized_blood_output.vcf.gz normalized_mom_output.vcf.gz

#bcftools annotate \
 #   -a "$cosmic" \
  #  -c CHROM,POS,INFO \
   # -o annotated_tumor_somatic_new38.vcf.gz \
   # -O z \
   # tumor_somatic_hg38.vcf.gz

#bcftools annotate \
 #   -a "$cosmic" \
  #  -c CHROM,POS,INFO \
   # -o annotated_shared_germline_new38.vcf.gz \
   # -O z \
   # blood_mom_somatic_hg38.vcf.gz





