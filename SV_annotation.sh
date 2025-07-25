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
bed="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RefSeq_Liftoff_v5.2.bed"
#nanomon_un="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv-new/tumor-all.nanomonsv.result.vcf"
nanomon="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv-new/tumor-all.nanomonsv.result.annotated.vcf"
nanomon_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv-new/tumor-all.nanomonsv.result.annrepeats.vcf"
sawfish="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic-new/sawfish-somatic.annotated.vcf"
sawfish_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish-somatic-new/sawfish-somatic.annrepeats.vcf"
savana="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana-new/savana.somatic.annotated.vcf"
savana_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-savana-new/savana.somatic.annrepeats.vcf"
severus="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-new/severus_somatic_breakpoints_double.annotated.vcf"
severus_annot="/lustre/scratch126/cellgen/behjati/lr26/PacBio-severus-new/severus_somatic_breakpoints_double.annrepeats.vcf"
repeats="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"

#bcftools sort "$nanomon_un" -o "$nanomon"

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

