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
somatic_tumor="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/isec_tumor_cosmic_somatic/0000.vcf.gz"
somatic_tumor_ann="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/annotated_tumor_somatic_cosmic_new38.vcf.gz"
### change to the tumor directory ###
cd "$dir"

### normalise all the vcf files using the reference genome (with decoys, same used for deepvariant calling) ###
bcftools norm -f "$reference" -Oz -o normalized_tumor_output.vcf.gz "$tumor"
bcftools norm -f "$reference" -Oz -o normalized_blood_output.vcf.gz "$blood"
bcftools norm -f "$reference" -Oz -o normalized_mom_output.vcf.gz "$mom"
bcftools norm -f "$reference" -Oz -o normalized_cosmic_hg38.vcf.gz "$cosmic"

### index all normalised vcf files ###
tabix -p vcf normalized_tumor_output.vcf.gz
tabix -p vcf normalized_blood_output.vcf.gz
tabix -p vcf normalized_mom_output.vcf.gz
tabix -p vcf normalized_cosmic_hg38.vcf.gz

### use bcftools to find all the variants that are only in the tumor file (somatic) ###
bcftools isec -p isec_tumor_cosmic_somatic -n=2 normalized_blood_output.vcf.gz normalized_mom_output.vcf.gz

### use bcftools annotate to annotate this file 
bcftools annotate \
    -a "$cosmic" \
    -c CHROM,POS,INFO \
    -o "$somatic_tumor_ann" \
    -O z \
    "$somatic_tumor"




