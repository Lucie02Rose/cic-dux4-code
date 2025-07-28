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
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic/tumor_somatic_annotated.vcf"
tumor_2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic/tumor_somatic_annotated_2.vcf"

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
### change to the output directory 
cd "$output_dir"

bcftools norm -m -any -Oz -o tumor_all_genotyped.norm.sv.vcf.gz "$tumor_all"
bcftools norm -m -any -Oz -o blood_genotyped.norm.sv.vcf.gz "$blood"

bcftools index tumor_all_genotyped.norm.sv.vcf.gz
bcftools index blood_genotyped.norm.sv.vcf.gz

# Unique to file1 (not in file2 or file3)
bcftools isec -p isec_tumor_only -C tumor_all_genotyped.norm.sv.vcf.gz blood_genotyped.norm.sv.vcf.gz


#cp "$nanomontum" . 

bcftools query -f '%CHROM\t%POS\t%INFO/END\n' sawfish_germline_blood_mom.vcf | awk '$3 < $2 {print $1"\t"$2}' > bad_sawfish_germline.txt
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' sawfish_tumor_somatic.vcf | awk '$3 < $2 {print $1"\t"$2}' > bad_sawfish_tumor.txt
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' tumor-somatic_nanomonsv.vcf | awk '$3 < $2 {print $1"\t"$2}' > bad_somatic_nanomonsv.txt

bcftools view -T ^bad_sawfish_germline.txt sawfish_germline_blood_mom.vcf -o valid_sawfish_germline_blood_mom.vcf -O v
bcftools view -T ^bad_sawfish_tumor.txt sawfish_tumor_somatic.vcf -o valid_sawfish_tumor_somatic.vcf -O v
bcftools view -T ^bad_somatic_nanomonsv.txt tumor_somatic_nanomonsv.vcf -o valid_tumor_somatic_nanomonsv.vcf -O v

bcftools view -T bad_sawfish_germline.txt sawfish_germline_blood_mom.vcf | bcftools annotate --remove INFO/END -o fixed_sawfish_germline_blood_mom.vcf -O v
bcftools view -T bad_sawfish_tumor.txt sawfish_tumor_somatic.vcf | bcftools annotate --remove INFO/END -o fixed_sawfish_tumor_somatic.vcf -O v
bcftools view -T bad_somatic_nanomonsv.txt tumor_somatic_nanomonsv.vcf | bcftools annotate --remove INFO/END -o fixed_tumor_somatic_nanomonsv.vcf -O v

bcftools concat -a valid_sawfish_germline_blood_mom.vcf fixed_sawfish_germline_blood_mom.vcf -o cleaned_sawfish_germline_blood_mom.vcf -O v
bcftools concat -a valid_sawfish_tumor_somatic.vcf fixed_sawfish_tumor_somatic.vcf -o cleaned_sawfish_tumor_somatic.vcf -O v
bcftools concat -a valid_tumor_somatic_nanomonsv.vcf fixed_tumor_somatic_nanomonsv.vcf -o cleaned_tumor_somatic_nanomonsv.vcf -O v

bcftools annotate -a "$bed" -c CHROM,FROM,TO,INFO/Gene -h <(echo '##INFO=<ID=Gene,Number=1,Type=String,Description="Overlapping gene">') -o annotated_sawfish_germline_blood_mom.vcf -O v cleaned_sawfish_germline_blood_mom.vcf 
bcftools annotate -a "$bed" -c CHROM,FROM,TO,INFO/Gene -h <(echo '##INFO=<ID=Gene,Number=1,Type=String,Description="Overlapping gene">') -o annotated_sawfish_tumor_somatic.vcf -O v cleaned_sawfish_tumor_somatic.vcf
bcftools annotate -a "$bed" -c CHROM,FROM,TO,INFO/Gene -h <(echo '##INFO=<ID=Gene,Number=1,Type=String,Description="Overlapping gene">') -o annotated_tumor_somatic_nanomonsv.vcf -O v cleaned_tumor_somatic_nanomonsv.vcf

