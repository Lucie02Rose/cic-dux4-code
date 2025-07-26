#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 300000
#BSUB -R 'span[hosts=1] select[mem>300000] rusage[mem=300000]'
#BSUB -q hugemem
#BSUB -J rna-t2t-alignment
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### Define directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"
output_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"
reference_index="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t/T2T_index"

# Move to the input directory
cd "$input_dir"

# 1️⃣ Generate STAR genome index (only if it does not exist)
if [ ! -f "$reference_index/Genome" ]; then
    echo "Generating STAR genome index..."
    STAR --runThreadN 16 \
         --runMode genomeGenerate \
         --genomeDir "$reference_index" \
         --genomeFastaFiles "$reference"
else
    echo "STAR genome index already exists. Skipping index generation."
fi

# 2️⃣ Function to run STAR alignment and BAM sorting
run_star_alignment() {
    sample=$1
    r1="$input_dir/${sample}_R1.fastq"
    r2="$input_dir/${sample}_R2.fastq"
    output_prefix="$output_dir/${sample}_"
    unsorted_bam="${output_prefix}Aligned.out.bam"
    sorted_bam="${output_prefix}Aligned.sortedByCoord.bam"

    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "Running STAR for paired-end reads: $sample"
        STAR --runThreadN 16 \
             --genomeDir "$reference_index" \
             --readFilesIn "$r1" "$r2" \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix "$output_prefix"

        # 3️⃣ Sort BAM file explicitly
        echo "Sorting BAM file for $sample..."
        samtools sort -@ 16 -o "$sorted_bam" "$unsorted_bam"

        # Remove unsorted BAM to save space
        rm "$unsorted_bam"
    else
        echo "ERROR: Missing FASTQ file(s) for sample $sample! Skipping..."
    fi
}

# 4️⃣ Find all paired-end samples and process in parallel (limit to 16 concurrent jobs)
for r1_file in "$input_dir"/*_R1.fastq; do
    sample_name=$(basename "$r1_file" | sed 's/_R1.fastq//')
    run_star_alignment "$sample_name" &

    # Limit concurrent jobs to 16
    if (( $(jobs | wc -l) >= 16 )); then
        wait -n  # Wait for at least one job to finish before launching a new one
    fi
done
wait  # Ensure all background jobs finish

echo "All STAR alignments and BAM sorting completed successfully!"


