#!/bin/bash
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-rna-our.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-rna-our.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J rna-t2t
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories to process
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
index="/lustre/scratch126/cellgen/behjati/lr26/T2T/salmon_T2T_index_te"
salmon="/lustre/scratch126/cellgen/behjati/lr26/RNA/Salmon"
#reference_index="/lustre/scratch126/cellgen/behjati/lr26/T2T/T2T_index"

# Define base directory where BAM files are stored
cram_file="/lustre/scratch126/cellgen/behjati/lr26/RNA/our_patient.cram"

# Create the output directory if it doesn't exist
cd "$output_dir"

sample_name="our_patient"

# Define output FASTQ file paths
r1_out="${output_dir}/${sample_name}_R1.fastq"
r2_out="${output_dir}/${sample_name}_R2.fastq"

# Generate FASTQ files from BAM using samtools
if [[ -f "$cram_file" ]]; then
    echo "Processing $sample_name..."
    samtools fastq -@ 16 -1 "$r1_out" -2 "$r2_out" "$cram_file"
    echo "Finished extracting FASTQ for $sample_name"
else
    echo "Warning: BAM file for $sample_name not found, skipping..."
fi

conda deactivate
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate salmon

salmon quant \
  -i "$index" \
  -l A \
  -1 "$r1_out" \
  -2 "$r2_out" \
  -p 16 \
  --validateMappings \
  -o "${salmon}/quant_te${sample_name}"

echo "Salmon quantification finished for $sample_name"

#conda deactivate
#source /software/cellgen/team274/lr26/miniforge3/bin/activate
#conda activate star_env

#sample_name="our_patient"
#r1="$r1_out"
#r2="$r2_out"
#output_prefix="$output_dir/${sample_name}_"
#unsorted_bam="${output_prefix}Aligned.out.bam"
#sorted_bam="${output_prefix}Aligned.sortedByCoord.bam"

# STAR alignment (only if FASTQ files exist)
#if [[ -f "$r1" && -f "$r2" ]]; then
   # echo "Running STAR for paired-end reads: $sample_name"
    
   # STAR --runThreadN 16 \
    #     --genomeDir "$reference_index" \
    #     --readFilesIn "$r1" "$r2" \
    #     --outSAMtype BAM Unsorted \
    #     --outFileNamePrefix "$output_prefix" \
    #     --chimSegmentMin 10 \
    #     --chimJunctionOverhangMin 20 \
    #     --chimOutType Junctions SeparateSAMold

   # echo "Sorting BAM file for $sample_name..."
   # samtools sort -@ 16 -o "$sorted_bam" "$unsorted_bam"

   # echo "Indexing sorted BAM file for $sample_name..."
   # samtools index "$sorted_bam"

   # rm "$unsorted_bam"

#else
   # echo "ERROR: Missing FASTQ file(s) for sample $sample_name! Skipping STAR..."
#fi
