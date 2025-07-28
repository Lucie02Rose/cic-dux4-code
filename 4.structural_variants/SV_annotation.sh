#!/bin/bash
### parameters for the LSF job ###
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-svfilt.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-svfilt.e
#BSUB -n 12
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J svfilt
#BSUB -G team274

### activate the conda envrionment with bcftools, samtools etc ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define all the directories - before running this script I ensure that I have all the output samples named ###
### how I want them in the correct directories - this means that I have moverd them to the respective directories here ###
### which is vital for this step to work ###
### this is also provided that sawfish filtering was carried out ###
### if there are any problems with these steps, they can be fixed by consulting the output and error files and using bcftools to sanitize variants ###
nanomon="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.sorted.vcf"
nanomon_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.annotated.vcf"
nanomon_repeats="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.annotated.repeats.vcf"
sawfish="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic-new/sawfish-somatic.annotated.vcf"
sawfish_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic-new/sawfish-somatic.annrepeats.vcf"
savana="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana-new/savana.somatic.annotated.vcf"
savana_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana-new/savana.somatic.annrepeats.vcf"
severus="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-new/severus_somatic_breakpoints_double.annotated.vcf"
severus_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-new/severus_somatic_breakpoints_double.annrepeats.vcf"
repeats="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
gff3_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
gff3="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"

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

### annotation steps for all the outputs ###
bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$nanomon_annot" \
  -O z \
  "$nanomon"

bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$savana_annot" \
  -O z \
  "$savana"

bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$severus_annot" \
  -O z \
  "$severus"

bcftools annotate \
  -a "$repeats" \
  -c CHROM,FROM,TO,INFO/repeats \
  -h <(echo '##INFO=<ID=repeats,Number=1,Type=String,Description="Repeats from T2T annotation">') \
  -o "$sawfish_annot" \
  -O z \
  "$sawfish"

