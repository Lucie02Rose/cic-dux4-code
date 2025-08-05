#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 32
#BSUB -M 200000
#BSUB -R 'span [hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q long
#BSUB -J salmon-cohort
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-salmon.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-salmon.e

### activate the conda bioinfo environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate salmon

### define all directoriess to process ###
### references ###
rna_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/GCF_009914755.1_T2T-CHM13v2.0_rna.fna.gz"
rna="/lustre/scratch126/cellgen/behjati/lr26/T2T/refseq_transcripts.fa"
### input, output index ###
input_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
index="/lustre/scratch126/cellgen/behjati/lr26/T2T/salmon_T2T_index_te"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA/Salmon"

### handling and indexing reference ###
gunzip -c "$rna_gz" > "$rna"
salmon index -t "$rna" -i "$index" --keepDuplicates

### making the output, changing to input ###
mkdir -p "$output_dir"
cd "$input_dir"

### processing fastq or fastq.gz sequentially making sure all are paired ###
for r1 in *_R1.fastq*; do
    r2="${r1/_R1.fastq/_R2.fastq}"
    ### extract sample name regardless of compression ###
    sample=$(basename "$r1" | sed -E 's/_R1\.fastq(.gz)?//')
    ### check if R2 is there and then run salmon quant ###
    if [[ -f "$r2" ]]; then
        salmon quant \
          -i "$index" \
          -l A \
          -1 "$r1" \
          -2 "$r2" \
          -p 32 \
          --validateMappings \
          -o "$output_dir/quant_te${sample}"
    else
        ### error handling ###
        echo "WARNING: No R2 file found for $r1. Skipping..."
    fi
done
