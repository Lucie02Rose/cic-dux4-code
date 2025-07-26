#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J pepper-blood
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/vepgermline.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/vepgermline.err

# Load correct compiler
export PATH=/usr/bin:$PATH
export PATH=$PATH:/nfs/users/nfs_l/lr26/ensembl-vep

# Define pathways
vep="/nfs/users/nfs_l/lr26/ensembl-vep"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/blood_common.vcf.gz"
dbsnp="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
gff="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gff3.gz"
clinvar="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0_ClinVar20220313.vcf.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01"

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
  --output_file "$output_dir/germline_annotated_clinvar_dbsnp.vcf" \
  --vcf \
  --force_overwrite

