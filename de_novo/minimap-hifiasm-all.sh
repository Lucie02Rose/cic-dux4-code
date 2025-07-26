#!/bin/bash
#BSUB -n 16
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q hugemem
#BSUB -J hifi
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/output_logs/hifi.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/error_logs/hifi.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.mmi"
samples=("PacBio-1B01-denovo")

# Loop over samples
for sample in "${samples[@]}"; do
    base_dir="/lustre/scratch126/casm/team274sb/lr26/"
    output_dir="$base_dir"
    tmp_dir="$output_dir/tmp"
    mkdir -p "$tmp_dir"
    export TMPDIR="$tmp_dir"

    # Find matching p_ctg FASTA files, excluding noseq and utg
    find "$base_dir" -maxdepth 1 -type f -name "${sample}_denovo_hifiasm.bp*.p_ctg.fasta" ! -name "*noseq*" ! -name "*utg*" | while read -r fasta; do
        # Extract filename without path
        filename=$(basename "$fasta")
        
        # Build output BAM name
        outname="${filename%.fasta}_vs_t2t"
        bam_out="${output_dir}/${outname}.bam"
   
        # Launch minimap2 in the background, using 16 cores
        echo "Processing $filename"
        (   
            minimap2 -t 16 -ax asm5 "$reference" "$fasta" | \
            samtools sort -@ 16 -o "$bam_out" - 
        
            samtools index "$bam_out"
        ) &
    done

    # Wait for all background jobs to finish before continuing
    wait
done



