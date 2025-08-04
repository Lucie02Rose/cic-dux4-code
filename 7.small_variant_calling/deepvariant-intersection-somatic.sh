#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J deepvariant
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-deepvar-intersec.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-deepvar-intersec.err

### activate the bioinfo conda environment with bcftools in it ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define t2t pathways ###
mom_vcf="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom/mom_output.vcf.gz"
blood_vcf="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood/blood_output.vcf.gz"
tumor_vcf="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor/tumor_output.vcf.gz"
mom_norm="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom/mom_output_norm.vcf.gz"
blood_norm="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood/blood_output_norm.vcf.gz"
tumor_norm="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor/tumor_output_norm.vcf.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor"
### define hg38 pathways ###
mom_vcf_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom-hg38/mom_output.vcf.gz"
blood_vcf_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood-hg38/blood_output.vcf.gz"
tumor_vcf_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/tumor_output.vcf.gz"
mom_norm_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom-hg38/mom_output_norm.vcf.gz"
blood_norm_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood-hg38/blood_output_norm.vcf.gz"
tumor_norm_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/tumor_output_norm.vcf.gz"
output_dir_hg38="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38"

### change to the tumor t2t reference output ###
cd "$output_dir"
### use bcftools to normalise and index all the vcf files ###
bcftools norm -m -any -Oz -o "$mom_norm" "$mom_vcf"
bcftools norm -m -any -Oz -o "$blood_norm" "$blood_vcf"
bcftools norm -m -any -Oz -o "$tumor_norm" "$tumor_vcf"
bcftools index "$mom_norm"
bcftools index "$blood_norm"
bcftools index "$tumor_norm"

### since the three have the same seq depth from the same platform I subtracted all that is present in the blood and mother ###
### what remains are the somatic tumor variants ###
bcftools isec -p isec1 -C "$tumor_norm" "$blood_norm" "$mom_norm"

### change to the tumor hg38 reference output ###
cd "$output_dir_hg38"
### use bcftools to normalise and index all the vcf files ###
bcftools norm -m -any -Oz -o "$mom_norm_hg38" "$mom_vcf_hg38"
bcftools norm -m -any -Oz -o "$blood_norm_hg38" "$blood_vcf_hg38"
bcftools norm -m -any -Oz -o "$tumor_norm_hg38" "$tumor_vcf_hg38"
bcftools index "$mom_norm_hg38"
bcftools index "$blood_norm_hg38"
bcftools index "$tumor_norm_hg38"

### since the three have the same seq depth from the same platform I subtracted all that is present in the blood and mother ###
### what remains are the somatic tumor variants ###
bcftools isec -p isec1 -C "$tumor_norm_hg38" "$blood_norm_hg38" "$mom_norm_hg38"
