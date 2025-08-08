#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 12
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J svannot
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-svannot.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-svannot.e

### activate the conda envrionment with bcftools, samtools etc ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define all the directories - before running this script I ensure that I have all the output samples named ###
### how I want them in the correct directories - this means that I have moverd them to the respective directories here ###
### which is vital for this step to work ###
### this is also provided that sawfish filtering and sorting of nanomonsv was carried out ###
### if there are any problems with these steps, they can be fixed by consulting the output and error files and using bcftools to sanitize variants ###
nanomon="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.sorted.vcf"
nanomon_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.annotated.vcf"
nanomon_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.annotated.repeats.vcf"
sawfish="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/isec_tumor_only/0000.vcf"
sawfish_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/sawfish-somatic.annotated.vcf"
sawfish_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/sawfish-somatic.annotated.repeats.vcf"
savana="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana/savana.somatic.vcf"
savana_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana/savana.somatic.annotated.vcf"
savana_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana/savana.somatic.annotated.repeats.vcf"
### note that for severus, it has firstly been converted from the tsv ouput to a vcf using the variant_callers_filtering.ipynb script ###
severus="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus/severus_somatic_breakpoints_double.vcf"
severus_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus/severus_somatic_breakpoints_double.annotated.vcf"
severus_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus/severus_somatic_breakpoints_double.annotated.repeats.vcf"
### sawfish germline also processed ###
sawfish_germline="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/blood_genotyped.norm.sort.sv.vcf.gz"
sawfish_germline_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/blood_genotyped.annotated.norm.sv.vcf.gz"
sawfish_germline_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/blood_genotyped.annotated.repeats.norm.sv.vcf.gz"
### files to use for annotation ###
repeats="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
gff3_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
gff3="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"

### convert the gff3 to the bed format for annotation, provided it is sorted by chromosome, extract what is needed ###
### and directly sort the bed file by chromosome ###
### repeat annotation from repeat masker does not really need annotation ###
echo "Converting gff3 file to a bed file "
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
echo "Indexing the bed file (optional)"
bgzip -c "$bed" > "${bed}.gz"
tabix -p bed "${bed}.gz"

### annotation steps for all the outputs ###
### nanomonsv ###
### bed file genic elements ###
echo "Annotating nanomonsv"
bcftools annotate \
  -a "$bed" \
  -c CHROM,FROM,TO,INFO/genes \
  -h <(echo '##INFO=<ID=genes,Number=1,Type=String,Description="Genes from T2T annotation">') \
  -o "$nanomon_annot" \
  -O z \
  "$nanomon"
### bed file repeat masker ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$nanomon_repeats" \
  -O z \
  "$nanomon_annot"
  
### savana ###
### bed file genic elements ###
echo "Annotating savana"
bcftools annotate \
  -a "$bed" \
  -c CHROM,FROM,TO,INFO/genes \
  -h <(echo '##INFO=<ID=genes,Number=1,Type=String,Description="Genes from T2T annotation">') \
  -o "$savana_annot" \
  -O z \
  "$savana"
### bed file repeat masker ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$savana_repeats" \
  -O z \
  "$savana_annot"
  
### severus ###
### bed file genic elements ###
echo "Annotating severus"
bcftools annotate \
  -a "$bed" \
  -c CHROM,FROM,TO,INFO/genes \
  -h <(echo '##INFO=<ID=genes,Number=1,Type=String,Description="Genes from T2T annotation">') \
  -o "$severus_annot" \
  -O z \
  "$severus"
### bed file repeat masker ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$severus_repeats" \
  -O z \
  "$severus_annot"
  
### sawfish ###
### bed file genic elements ###
echo "Annotating sawfish"
bcftools annotate \
  -a "$bed" \
  -c CHROM,FROM,TO,INFO/genes \
  -h <(echo '##INFO=<ID=genes,Number=1,Type=String,Description="Genes from T2T annotation">') \
  -o "$sawfish_annot" \
  -O z \
  "$sawfish"
### bed file repeat masker ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$sawfish_repeats" \
  -O z \
  "$sawfish_annot"

### sawfish germline ###
### bed file genic elements ###
echo "Annotating sawfish"
bcftools annotate \
  -a "$bed" \
  -c CHROM,FROM,TO,INFO/genes \
  -h <(echo '##INFO=<ID=genes,Number=1,Type=String,Description="Genes from T2T annotation">') \
  -o "$sawfish_germline_annot" \
  -O z \
  "$sawfish_germline"
### bed file repeat masker ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$sawfish_germline_repeats" \
  -O z \
  "$sawfish_germline_annot"

