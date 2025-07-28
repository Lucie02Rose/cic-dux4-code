#!/bin/bash
### details for LSF job ###
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-hifi-new-hg38
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J_pbmm2-hifi-new-hg38.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J_pbmm2-hifi-new-hg38.err

### activate my base conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

#### directories including temporary for sorting ###
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned-hg38"
tmp_dir="$output_dir/tmp"
reference_fasta="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.fa"
reference_index="/lustre/scratch126/cellgen/behjati/lr26/hg38/hg38.mmi"

### check if the index or reference exists, if not then make it ###
### echo are used as messages written to the -o -e files ###
if [ ! -f "$reference_index" ]; then
    echo "Index not found at $reference_index. Creating pbmm2 index..."
    pbmm2 index "$reference_fasta" "$reference_index"
    echo "Indexing completed."
else
    echo "Reference index already exists at $reference_index."
fi

### create the output and temporary directories if not alreads there ###
mkdir -p "$output_dir" "$tmp_dir"
### pbmm2 needs to export the temporary for sorting ###
export TMPDIR="$tmp_dir"

### a for loop for the bam files in the input to ###
for input_bam in "$input_dir"/*.bam; do
    ### extract the base name without the .bam for all samples ###
    base_name=$(basename "$input_bam" .bam)
    ### define what will the output look like and give it the _pbmm2.bam extension ###
    output_bam="$output_dir/${base_name}_pbmm2.bam"
    ### print statement for the files being aligned ###
    echo "Aligning $input_bam to $reference..."
    ### pbmm2 align with the preset HIFI for high fidelity reads, sort for direct sorting, j for multithreaded, unmapped for mapping 
    pbmm2 align "$reference_index" "$input_bam" "$output_bam" --preset HIFI --sort -j 16 --unmapped
    ### echo for having finished ###
    echo "pbmm2 alignment completed: $output_bam"
done
