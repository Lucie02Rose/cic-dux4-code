#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J pbmm2-giab
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-maf.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-maf.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate wgatools

cd /lustre/scratch126/cellgen/behjati/lr26/PacBio-last/ 

# Directories
tumor_1B01="/lustre/scratch126/cellgen/behjati/lr26/PacBio-last/tumor_1B01_hifi_reads.maf"
tumor_extract="/lustre/scratch126/cellgen/behjati/lr26/PacBio-last/tumor_1B01_hifi_reads-extracted.maf"

grep -A 3 "^s chr4" "$tumor_1B01" > filtered_chr4.maf
grep -A 3 "^s chr10" "$tumor_1B01" > filtered_chr10.maf
cat filtered_chr4.maf filtered_chr10.maf > tumor_chr4_chr10_only.maf

wgatools maf-index tumor_chr4_chr10_only.maf

wgatools me tumor_chr4_chr10_only.maf \
  -r chr4:193430000-193570000,chr10:134610000-134750000 \
  > "tumor_extract"

