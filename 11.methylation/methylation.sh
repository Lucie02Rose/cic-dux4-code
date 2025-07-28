#!/bin/bash
### job specifications for the LSF ###
#BSUB -n 16
#BSUB -M 6000
#BSUB -R 'span[hosts=1] select[mem>6000] rusage[mem=6000]'
#BSUB -q long
#BSUB -J methylation
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-methylation.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-methylation.err

### activate conda environment for methylation ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate methylation

### directories including the reference - which can be the mmi index here ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.mmi"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-methylation"
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"

### create the methylation output directory ###
mkdir -p "$output_dir"

### message for the output ###
echo "Running methylation on all BAM files in $input_dir..."

### finding all files in the output that end with bam, retain their sample name
for input_bam in "$input_dir"/*.bam; do
    sample_name=$(basename "$input_bam" .bam)
    ### message - which sample is processed ###
    echo "Processing: $input_bam"
    ### running the methylation with 16 threads in the for loop one sample after the other ###
    aligned_bam_to_cpg_scores --bam "$input_bam" --output-prefix  "$output_dir/${sample_name}-methylation" --threads 16
    ### display message when file is done ###
    echo "Completed: $input_bam"
done
### display message when all is done ###
echo "All methylation analyses completed."


