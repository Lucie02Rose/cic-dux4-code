#!/bin/bash
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-salmon.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-salmon.e
#BSUB -n 32
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J pbmm2indexing
#BSUB -G team274

source /software/cellgen/team274/lr26/miniforge3/bin/activate
#conda activate bioinfo

### Directories to process
rna_gz="/lustre/scratch126/cellgen/behjati/lr26/T2T/GCF_009914755.1_T2T-CHM13v2.0_rna.fna.gz"
rna="/lustre/scratch126/cellgen/behjati/lr26/T2T/refseq_transcripts.fa"
repeats="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
repeats_fa="/lustre/scratch126/cellgen/behjati/lr26/T2T/repeatmasker_te.fa"
combined_transcripts="/lustre/scratch126/cellgen/behjati/lr26/T2T/combined_transcripts.fa"
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
index="/lustre/scratch126/cellgen/behjati/lr26/T2T/salmon_T2T_index_te"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA/Salmon"

#bedtools getfasta -fi "$reference" -bed "$repeats" -s -name > "$repeats_fa"
#gunzip -c "$rna_gz" > "$rna"
#cat "$rna" "$repeats_fa" > "$combined_transcripts"

#conda deactivate
conda activate salmon

# Indexing the T2T reference genome (if index files don't exist)
salmon index -t "$rna" -i "$index" --keepDuplicates

#mkdir -p "$output_dir"
cd "$input_dir"

for r1 in *FO_R1.fastq*; do
    # Derive the corresponding R2 file
    r2="${r1/_R1.fastq/_R2.fastq}"
    
    # Extract sample name (handles .fastq and .fastq.gz)
    sample=$(basename "$r1" | sed -E 's/_R1\.fastq(.gz)?//')
    
    # Make sure R2 file exists
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
        echo "WARNING: No R2 file found for $r1. Skipping..."
    fi
done
