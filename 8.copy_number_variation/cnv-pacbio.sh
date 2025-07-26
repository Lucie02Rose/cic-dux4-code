#!/bin/bash
#BSUB -n 16
#BSUB -M 150000
#BSUB -R 'span[hosts=1] select[mem>150000] rusage[mem=150000]'
#BSUB -q basement
#BSUB -J cnv-pacbio
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-cnv-pacbio.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-cnv-pacbio.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
ref="/lustre/scratch126/cellgen/behjati/lr26/T2T/"
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
dbsnp="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
dbsnp_sub="/lustre/scratch126/cellgen/behjati/lr26/T2T/t2t_chr1_10_dbSNPv155.vcf.gz"
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
mom="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/mom_1B02_hifi_reads_pbmm2.bam"
patient="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/1C01_blood_hifi_reads_pbmm2.bam"
dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-hificnv-tumor-all"

cd "$ref"

#sample only dbsnps on chr1 and 10 and index the vcf
bcftools view -r chr1,chr10 "$dbsnp" -Oz -o "$dbsnp_sub"
bcftools index "$dbsnp_sub"

cd "$dir"

bcftools mpileup -f "$reference" -a FORMAT/AD -l "$dbsnp_sub" -o tumor.bcf "$tumor"
bcftools view tumor.bcf -Ov -o tumor.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' tumor.vcf | \
awk -F'\t' '{
    split($5, ad, ","); 
    total=0; 
    for (i in ad) total+=ad[i]; 
    if (total > 0) {
        printf "%s", $0;
        for (i in ad) printf "\t%.3f", ad[i] / total;
        print "";
    } else {
        print $0, "0";
    }
}' > tumor_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -l "$dbsnp_sub" -o blood.bcf "$blood"
bcftools view blood.bcf -Ov -o blood.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' blood.vcf | \
awk -F'\t' '{
    split($5, ad, ",");
    total=0;
    for (i in ad) total+=ad[i];
    if (total > 0) {
        printf "%s", $0;
        for (i in ad) printf "\t%.3f", ad[i] / total;
        print "";
    } else {
        print $0, "0";
    }
}' > blood_fin.vcf
