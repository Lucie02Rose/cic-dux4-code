#!/bin/bash
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-rna-star.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-rna-star.e
#BSUB -n 64
#BSUB -M 600000
#BSUB -R 'span[hosts=1] select[mem>600000] rusage[mem=600000]'
#BSUB -q hugemem
#BSUB -J rna-t2t-alignment
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### Define directories
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/RNA/STAR-T2T"
reference_index="/lustre/scratch126/cellgen/behjati/lr26/T2T/T2T_index"

mkdir -p "$output_dir"
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

run_star_alignment() {
    sample=$1

    # Try both fastq and fastq.gz extensions
    r1="$input_dir/${sample}_R1.fastq"
    r2="$input_dir/${sample}_R2.fastq"
    if [[ ! -f "$r1" ]]; then
        r1="$input_dir/${sample}_R1.fastq.gz"
    fi
    if [[ ! -f "$r2" ]]; then
        r2="$input_dir/${sample}_R2.fastq.gz"
    fi

    output_prefix="$output_dir/${sample}_"
    unsorted_bam="${output_prefix}Aligned.out.bam"
    sorted_bam="${output_prefix}Aligned.sortedByCoord.bam"

    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "Running STAR for paired-end reads: $sample"
        STAR --runThreadN 16 \
             --genomeDir "$reference_index" \
             --readFilesIn "$r1" "$r2" \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix "$output_prefix" \
             --chimSegmentMin 10 \
             --chimJunctionOverhangMin 20 \
             --chimOutType Junctions SeparateSAMold

        echo "Sorting BAM file for $sample..."
        samtools sort -@ 16 -o "$sorted_bam" "$unsorted_bam"

        echo "Indexing sorted BAM file for $sample..."
        samtools index "$sorted_bam"

        rm "$unsorted_bam"
    else
        echo "ERROR: Missing FASTQ file(s) for sample $sample! Skipping..."
    fi
}


# Find all paired-end samples and process in parallel (limit to 16 concurrent jobs)
for r1_file in "$input_dir"/*_R1.fastq*; do
    # Determine the corresponding R2 filename
    r2_file="${r1_file/_R1.fastq/_R2.fastq}"
    r2_file="${r2_file/.gz/.gz}"  # ensures consistency if it's already gzipped

    if [[ -f "$r2_file" ]]; then
        # Extract sample name (removes _R1.fastq or _R1.fastq.gz)
        sample_name=$(basename "$r1_file" | sed -E 's/_R1\.fastq(.gz)?//')
        run_star_alignment "$sample_name" &
    else
        echo "WARNING: No R2 file found for $(basename "$r1_file"). Skipping..."
    fi

    # Limit concurrent jobs to 16
    if (( $(jobs | wc -l) >= 16 )); then
        wait -n
    fi
done
wait  # Ensure all background jobs finish


echo "All STAR alignments, BAM sorting, and indexing completed successfully!"
