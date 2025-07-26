#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J pepper-blood
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/vepblood.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/vepblood.err

export PATH=/usr/bin:$PATH
export PATH=$PATH:/nfs/users/nfs_l/lr26/ensembl-vep

# Define pathways
vep="/nfs/users/nfs_l/lr26/ensembl-vep"
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
input_vcf="/lustre/scratch126/casm/team274sb/lr26/pepper-blood38/blood_output.vcf.gz"
dbsnp="/lustre/scratch126/casm/team274sb/lr26/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
cosmic="/lustre/scratch126/casm/team274sb/lr26/hg38/cleaned_cosmic_hg38_nochr.vcf.gz"
gff="/lustre/scratch126/casm/team274sb/lr26/hg38/gencode.sorted.gff3.gz"
clinvar="/lustre/scratch126/casm/team274sb/lr26/hg38/clinvar.vcf.gz"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pepper-blood38"

# Change to pepper
cd "$output_dir"

# Run vep with custom ClinVar and dbSNP liftover annotations
vep -i "$input_vcf" \
    -o "$output_dir/blood_vep_annotated_cosmic_clinvar_dbsnp.vcf" \
    --fasta "$reference" \
    --gff "$gff" \
    --species homo_sapiens \
    --assembly GRCh38 \
    --vcf \
    --custom "$dbsnp",dbSNP,vcf,exact,0,ID \
    --custom "$clinvar",ClinVar,vcf,overlap,0,CLNSIG \
    --custom "$cosmic",cosmic,vcf,exact,0,COSMIC_ID \
    --force_overwrite
