#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J star-our
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-star-our.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-star-our.e

### activate conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process including the two downloaded and renamed cram files ###
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
reference_index="/lustre/scratch126/cellgen/behjati/lr26/T2T/T2T_index"
cram_FO="/lustre/scratch126/cellgen/behjati/lr26/RNA/our_patient_FO.cram"
cram_FT="/lustre/scratch126/cellgen/behjati/lr26/RNA/our_patient_FT.cram"

### change to the output directory ###
cd "$output_dir"
### define sample names and output fastq paths ###
FO_sample="our_patient_FO"
FT_sample="our_patient_FT"
r1_FO="${output_dir}/${FO_sample}_R1.fastq"
r2_FO="${output_dir}/${FO_sample}_R2.fastq"
r1_FT="${output_dir}/${FT_sample}_R1.fastq"
r2_FT="${output_dir}/${FT_sample}_R2.fastq"

### generate fastq files from the cram files ###
if [[ -f "$cram_FO" ]]; then
    echo "Processing $FO_sample..."
    samtools fastq -@ 16 -1 "$r1_FO" -2 "$r2_FO" "$cram_FO"
    echo "Finished extracting FASTQ for $FO_sample"
else
    echo "Warning: BAM file for $FO_sample not found, skipping..."
fi
### repeat for FT sample ###
if [[ -f "$cram_FT" ]]; then
    echo "Processing $FT_sample..."
    samtools fastq -@ 16 -1 "$r1_FT" -2 "$r2_FT" "$cram_FT"
    echo "Finished extracting FASTQ for $FT_sample"
else
    echo "Warning: BAM file for $FT_sample not found, skipping..."
fi

### now run star after activating the conda for star ###
conda deactivate
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### define outputs for both samples ###
output_prefix_FO="$output_dir/${FO_sample}_"
unsorted_bam_FO="${output_prefix_FO}Aligned.out.bam"
sorted_bam_FO="${output_prefix_FO}Aligned.sortedByCoord.bam"

output_prefix_FT="$output_dir/${FT_sample}_"
unsorted_bam_FT="${output_prefix_FT}Aligned.out.bam"
sorted_bam_FT="${output_prefix_FT}Aligned.sortedByCoord.bam"

### run star on the fo sample if paired end fastq are present ###
if [[ -f "$r1_FO" && -f "$r2_FO" ]]; then
    echo "Running STAR for paired-end reads: $FO_sample"
    ### star settings ###
    STAR --runThreadN 16 \
         --genomeDir "$reference_index" \
         --readFilesIn "$r1_FO" "$r2_FO" \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix "$output_prefix_FO" \
         --chimSegmentMin 10 \
         --chimJunctionOverhangMin 20 \
         --chimOutType Junctions SeparateSAMold
    ### sorting and indexing ###
    echo "Sorting BAM file for $FO_sample..."
    samtools sort -@ 16 -o "$sorted_bam_FO" "$unsorted_bam_FO"
    echo "Indexing sorted BAM file for $FO_sample..."
    samtools index "$sorted_bam_FO"
    ### remove unsorted ###
    rm "$unsorted_bam_FO"
else
    ### error handling ###
    echo "ERROR: Missing FASTQ file(s) for sample $FO_sample! Skipping STAR..."
fi

### run star on the fo sample if paired end fastq are present ###
if [[ -f "$r1_FT" && -f "$r2_FT" ]]; then
    echo "Running STAR for paired-end reads: $FT_sample"
    ### star settings ###
    STAR --runThreadN 16 \
         --genomeDir "$reference_index" \
         --readFilesIn "$r1_FT" "$r2_FT" \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix "$output_prefix_FT" \
         --chimSegmentMin 10 \
         --chimJunctionOverhangMin 20 \
         --chimOutType Junctions SeparateSAMold
    ### sorting and indexing ###
    echo "Sorting BAM file for $FT_sample..."
    samtools sort -@ 16 -o "$sorted_bam_FT" "$unsorted_bam_FT"
    echo "Indexing sorted BAM file for $FT_sample..."
    samtools index "$sorted_bam_FT"
    ### remove unsorted ###
    rm "$unsorted_bam_FT"
else
    ### error handling ###
    echo "ERROR: Missing FASTQ file(s) for sample $FT_sample! Skipping STAR..."
fi
