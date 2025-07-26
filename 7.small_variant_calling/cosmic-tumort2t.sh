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
cosmict2t="/lustre/scratch126/casm/team274sb/lr26/T2T/normalized_cosmic_T2T.vcf.gz"
referencet2t="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
reference38="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
tumor38="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor38/tumor_output.vcf.gz"
cosmic38="/lustre/scratch126/casm/team274sb/lr26/hg38/cleaned_cosmic_hg38.vcf.gz"
tumort2t="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01/tumor_output.vcf.gz"
dir="/lustre/scratch126/casm/team274sb/lr26/pepper-tumor1B01"

cd "$dir"

bcftools norm -f "$reference38" -Oz -o normalized_tumor_outputhg38.vcf.gz "$tumor38"
bcftools norm -f "$referencet2t" -Oz -o normalized_tumor_outputt2t.vcf.gz "$tumort2t"

tabix -p vcf normalized_tumor_outputhg38.vcf.gz
tabix -p vcf normalized_tumor_outputt2t.vcf.gz

bcftools annotate \
    -a "$cosmict2t" \
    -c CHROM,POS,INFO \
    -o annotated_normalized_tumor_outputt2t.vcf.gz \
    -O z \
    normalized_tumor_outputt2t.vcf.gz

bcftools annotate \
    -a "$cosmic38" \
    -c CHROM,POS,INFO \
    -o annotated_normalized_tumor_outputhg38.vcf.gz \
    -O z \
    normalized_tumor_outputhg38.vcf.gz





