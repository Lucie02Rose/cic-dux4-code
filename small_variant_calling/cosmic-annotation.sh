#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J cosmic
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/annot.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/annot.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
cosmic="/lustre/scratch126/casm/team274sb/lr26/T2T/Cosmic_MutantCensus_sorted_T2T.vcf.gz"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
tumor_somatic="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/isec1/tumor_somatic.vcf.gz"
common_germline="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/isec_2_3_common/blood_common.vcf.gz"
dir="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01"

cd "$dir"

bcftools norm -f "$reference" -Oz -o normalized_tumor_somatic.vcf.gz "$tumor_somatic"
bcftools norm -f "$reference" -Oz -o normalized_shared_germline.vcf.gz "$common_germline"
bcftools norm -f "$reference" -Oz -o normalized_cosmic_T2T.vcf.gz "$cosmic"

tabix -p vcf normalized_tumor_somatic.vcf.gz
tabix -p vcf normalized_shared_germline.vcf.gz
tabix -p vcf normalized_cosmic_T2T.vcf.gz

bcftools annotate \
    -a normalized_cosmic_T2T.vcf.gz \
    -c CHROM,POS,INFO \
    -o normalized_annotated_tumor_somatic_new.vcf.gz \
    -O z \
    normalized_tumor_somatic.vcf.gz

bcftools annotate \
    -a normalized_cosmic_T2T.vcf.gz \
    -c CHROM,POS,INFO \
    -o normalized_annotated_shared_germline_new.vcf.gz \
    -O z \
    normalized_shared_germline.vcf.gz





