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
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish"
tumor_all="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/tumor_all_4_hifi_reads_pbmm2_joint_call/genotyped.sv.vcf.gz"
blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/blood_1C01_hifi_reads_pbmm2_joint_call/genotyped.sv.vcf.gz"
gff3_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
gff3="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/tumor_somatic_annotated.vcf"

### convert the gff3 to the bed format for annotation, provided it is sorted by chromosome, extract what is needed ###
### and directly sort the bed file by chromosome ###
zcat "$gff3_gz" | awk -F'\t' '
$3 ~ /^(gene|transcript|exon|CDS|five_prime_UTR|three_prime_UTR)$/ {
    id = "."  # default
    if (match($9, /gene_name=([^;]+)/, a)) {
        id = a[1]
    } else if (match($9, /ID=([^;]+)/, a)) {
        id = a[1]
    }
    print $1, $4 - 1, $5, id, ".", $7
}' OFS='\t' | sort -k1,1 -k2,2n > "$bed"
### zip and index the bed file to use for annotation ###
bgzip -c "$bed" > "${bed}.gz"
tabix -p bed "${bed}.gz"
### change to the output directory ###
cd "$output_dir"
### normalise the sawfish output ###
bcftools norm -m -any -Oz -o tumor_all_genotyped.norm.sv.vcf.gz "$tumor_all"
bcftools norm -m -any -Oz -o blood_genotyped.norm.sv.vcf.gz "$blood"
### sort the sawfish output ###
bcftools sort tumor_all_genotyped.norm.sv.vcf.gz -Oz -o tumor_all_genotyped.norm.sort.sv.vcf.gz
bcftools sort blood_genotyped.norm.sv.vcf.gz -Oz -o blood_genotyped.norm.sort.sv.vcf.gz
### index the sawfish output ###
bcftools index tumor_all_genotyped.norm.sort.sv.vcf.gz
bcftools index blood_genotyped.norm.sort.sv.vcf.gz
### filter what is unique to the tumor (keep somatic variants ###
bcftools isec -p isec_tumor_only -C tumor_all_genotyped.norm.sort.sv.vcf.gz blood_genotyped.norm.sort.sv.vcf.gz
### normalisation, sorting and indexing should fix any downstream issues ###
### however if there are errors e.g. malformed starts and end positions (actually might be an inversion) ###
### try querying the ones where the start iw bigger than the end, e.g. awk '$3 < $2 {print $1"\t"$2}' ###
### label those as bad variants in a txt file and then using view exclude all the ones from the file cleaning the file ####
