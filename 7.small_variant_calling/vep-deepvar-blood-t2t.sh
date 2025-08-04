#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J vep-blood
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-vep-blood.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-vep-blood.err

### load the correct compiler to work with vep ###
export PATH=/usr/bin:$PATH
export PATH=$PATH:/nfs/users/nfs_l/lr26/ensembl-vep

### input, vep and reference ###
vep="/nfs/users/nfs_l/lr26/ensembl-vep"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_vcf="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood/blood_output_norm.vcf.gz"
### files to annotate with ###
dbsnp="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
gff="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gff3.gz"
clinvar="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_ClinVar20220313.vcf.gz"
### output ###
output_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-blood"

### change to the output directory ###
cd "$output_dir"

### run vep with custom ClinVar and dbSNP liftover annotations ###
vep \
  --i "$input_vcf" \
  --fasta "$reference" \
  --gff "$gff" \
  --species homo_sapiens \
  --assembly T2T-CHM13v2.0 \
  --custom "$clinvar",ClinVar,vcf,overlap,0,CLNSIG \
  --custom "$dbsnp",dbSNP,vcf,exact,0,ID \
  --o "$output_dir/blood_vep_annotated_with_clinvar_and_dbsnp_vcf_new.vcf" \
  --vcf \
  --force_overwrite

