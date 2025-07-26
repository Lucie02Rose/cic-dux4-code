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
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-filtering"
tumor_all="/lustre/scratch126/cellgen/behjati/lr26/mor-all/tumor-all_pbmm2-farm22_joint_call/genotyped.sv.vcf.gz"
#mom="/lustre/scratch126/cellgen/behjati/lr26/sawfish/1B02genotyped.sv.vcf.gz"
blood="/lustre/scratch126/cellgen/behjati/lr26/sawfish/1C01genotyped.sv.vcf.gz"
gff3_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
gff3="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"

# decompress
gunzip -c "$gff3_gz" > "$gff3"

# convert GFF3 to BED
gffread "$gff3" -T -o - | awk 'BEGIN{OFS="\t"} $3=="transcript" {print $1, $4-1, $5, $10, ".", $7}' > "$bed"

#mkdir -p "$output_dir"
cd "$output_dir"

#bcftools norm -m -any -Oz -o tumor_all_genotyped.norm.sv.vcf.gz "$tumor_all"
#bcftools norm -m -any -Oz -o blood_genotyped.norm.sv.vcf.gz "$mom"
#bcftools norm -m -any -Oz -o mom_genotyped.norm.sv.vcf.gz "$blood"

#bcftools index tumor_all_genotyped.norm.sv.vcf.gz
#bcftools index blood_genotyped.norm.sv.vcf.gz
#bcftools index mom_genotyped.norm.sv.vcf.gz

# Unique to file1 (not in file2 or file3)
#bcftools isec -p isec_tumor_only -C tumor_all_genotyped.norm.sv.vcf.gz blood_genotyped.norm.sv.vcf.gz mom_genotyped.norm.sv.vcf.gz

#bcftools isec -p isec_germline_blood_mom -n=2 blood_genotyped.norm.sv.vcf.gz mom_genotyped.norm.sv.vcf.gz

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

