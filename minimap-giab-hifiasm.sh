#!/bin/bash
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q yesterday
#BSUB -J hifi
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-hifi-giab.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-hifi-giab.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"
samples=("PacBio-denovo")

# Loop over samples
for sample in "${samples[@]}"; do
    base_dir="/lustre/scratch126/cellgen/behjati/lr26/"
    output_dir="$base_dir"
    tmp_dir="$output_dir/tmp"
    mkdir -p "$tmp_dir"
    export TMPDIR="$tmp_dir"

    # Find matching p_ctg FASTA files, excluding noseq and utg
    find "$base_dir" -maxdepth 1 -type f -name "${sample}_new.bp*.p_ctg.fasta" ! -name "*noseq*" ! -name "*utg*" | while read -r fasta; do
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



