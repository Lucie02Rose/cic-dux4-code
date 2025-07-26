#!/bin/bash
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-svfilt.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-svfilt.e
#BSUB -n 12
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J svfilt
#BSUB -G team274

source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories to process
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic-new"
tumor_all="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-tumor-all/tumor_all_4_hifi_reads_pbmm2_joint_call/genotyped.sv.vcf.gz"
blood="/lustre/scratch126/cellgen/behjati/lr26/genotyped.sv.vcf.gz"
gff3_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
gff3="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic/tumor_somatic_annotated.vcf"
tumor_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic/tumor_somatic_annotated_2.vcf"

# decompress
#gunzip -c "$gff3_gz" > "$gff3"

# convert GFF3 to BED
#zcat "$gff3_gz" | awk -F'\t' '
#$3 ~ /^(gene|transcript|exon|CDS|five_prime_UTR|three_prime_UTR)$/ {
 #   id = "."  # default
  #  if (match($9, /gene_name=([^;]+)/, a)) {
   #     id = a[1]
   # } else if (match($9, /ID=([^;]+)/, a)) {
   #     id = a[1]
   # }
   # print $1, $4 - 1, $5, id, ".", $7
#}' OFS='\t' > "$bed"

#bgzip -c "$bed" > "${bed}.gz"
#tabix -p bed "${bed}.gz"

mkdir -p "$output_dir"
cd "$output_dir"

bcftools norm -m -any -Oz -o tumor_all_genotyped.norm.sv.vcf.gz "$tumor_all"
bcftools norm -m -any -Oz -o blood_genotyped.norm.sv.vcf.gz "$blood"

bcftools index tumor_all_genotyped.norm.sv.vcf.gz
bcftools index blood_genotyped.norm.sv.vcf.gz

# Unique to file1 (not in file2 or file3)
bcftools isec -p isec_tumor_only -C tumor_all_genotyped.norm.sv.vcf.gz blood_genotyped.norm.sv.vcf.gz


