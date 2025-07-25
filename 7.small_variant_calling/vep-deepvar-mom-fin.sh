#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J pepper-mom
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

export PATH=/usr/bin:$PATH
export PATH=$PATH:/nfs/users/nfs_l/lr26/ensembl-vep

# Define pathways
vep="/nfs/users/nfs_l/lr26/ensembl-vep"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-mom/mom_output.vcf.gz"
dbsnp="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
gff="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gff3.gz"
clinvar="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_ClinVar20220313.vcf.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pepper-mom"

# Change to pepper
cd "$output_dir"

# Run vep with custom ClinVar and dbSNP liftover annotations
vep \
  --input_file "$input_vcf" \
  --fasta "$reference" \
  --gff "$gff" \
  --species homo_sapiens \
  --assembly T2T-CHM13v2.0 \
  --custom "$clinvar",ClinVar,vcf,overlap,0,CLNSIG \
  --custom "$dbsnp",dbSNP,vcf,exact,0,ID \
  --output_file "$output_dir/mom_vep_annotated_with_clinvar_and_dbsnp_vcf_new.vcf" \
  --vcf \
  --force_overwrite

