#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J vep-tum-som
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-vep-tum-som-hg38.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-vep-tum-som-hg38.err

### load the correct compiler to work with vep ###
export PATH=/usr/bin:$PATH
export PATH=$PATH:/nfs/users/nfs_l/lr26/ensembl-vep

### input, vep and reference ###
vep="/nfs/users/nfs_l/lr26/ensembl-vep"
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/hg38.fa"
input_vcf="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38/isec1/0000.vcf.gz"
### files to annotate with ###
dbsnp="/lustre/scratch126/casm/team274sb/lr26/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
gff="/lustre/scratch126/casm/team274sb/lr26/hg38/gencode.sorted.gff3.gz"
clinvar="/lustre/scratch126/casm/team274sb/lr26/hg38/clinvar.vcf.gz"
### output ###
output_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-tumor-hg38"

### change to the output directory ###
cd "$output_dir"

### run vep with custom ClinVar and dbSNP liftover annotations ###
vep \
  --i "$input_vcf" \
  --fasta "$reference" \
  --gff "$gff" \
  --species homo_sapiens \
  --assembly GRCh38 \
  --custom "$clinvar",ClinVar,vcf,overlap,0,CLNSIG \
  --custom "$dbsnp",dbSNP,vcf,exact,0,ID \
  --o "$output_dir/tumor_vep_annotated_with_clinvar_and_dbsnp_vcf_new_somatic_hg38.vcf.gz" \
  --vcf \
  --force_overwrite
