#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J cosmic
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-cosmic.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-cosmic.err

### activate the bioinfo conda environment with bcftools ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define the directories (COSMIC vcf is from the database) ###
cosmic="/lustre/scratch126/casm/team274sb/lr26/hg38/cleaned_cosmic_hg38.vcf.gz"
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/hg38.fa"
mom="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom-hg38/mom_output.vcf.gz"
blood="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood-hg38/blood_output.vcf.gz"
tumor="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/tumor_output.vcf.gz"
dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38"

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

bcftools annotate \
    -a "$cosmic" \
    -c CHROM,POS,INFO \
    -o annotated_tumor_somatic_new38.vcf.gz \
    -O z \
    tumor_somatic_hg38.vcf.gz

bcftools annotate \
    -a "$cosmic" \
    -c CHROM,POS,INFO \
    -o annotated_shared_germline_new38.vcf.gz \
    -O z \
    blood_mom_somatic_hg38.vcf.gz





