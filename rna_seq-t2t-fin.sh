#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q hugemem
#BSUB -J rna-t2t-alignemnt
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### Directories to process
reference="/nfs/users/nfs_l/lr26/PacBio/T2T/chm13v2.0.fa"
output_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"
reference_index="lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t/T2T_index"

# Move to the output directory
cd "$output_dir"

# 1️⃣ Generate STAR genome index (only if it does not exist)
if [ ! -f "$output_dir/T2T_index/Genome" ]; then
    echo "Generating STAR genome index..."
    STAR --runThreadN 16 \
         --runMode genomeGenerate \
         --genomeDir "$output_dir/T2T_index" \
         --genomeFastaFiles "$reference"
else
    echo "STAR genome index already exists. Skipping index generation."
fi

# 2️⃣ Run STAR alignment for paired-end reads
if [ -f "reads_R1.fastq" ] && [ -f "reads_R2.fastq" ]; then
    echo "Running STAR for paired-end reads..."
    STAR --runThreadN 16 \
         --genomeDir "$output_dir/T2T_index" \
         --readFilesIn reads_R1.fastq reads_R2.fastq \
         --limitBAMsortRAM 100000000000 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMheaderHD @HD VN:1.6 SO:coordinate \
         --outFileNamePrefix "$output_dir/star_paired_"
else
    echo "ERROR: Paired-end FASTQ files not found!"
    exit 1
fi

# 3️⃣ Run STAR alignment for singleton reads (if file exists)
if [ -f "singletons.fastq" ]; then
    echo "Running STAR for singleton reads..."
    STAR --runThreadN 16 \
         --genomeDir "$output_dir/T2T_index" \
         --readFilesIn singletons.fastq \
         --limitBAMsortRAM 160000000000 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMheaderHD @HD VN:1.6 SO:coordinate \
	 --outFileNamePrefix "$output_dir/star_single_"
else
    echo "WARNING: Singleton FASTQ file not found! Skipping..."
fi

# 4️⃣ Merge BAM files if both exist
if [ -f "$output_dir/star_paired_Aligned.sortedByCoord.out.bam" ] && [ -f "$output_dir/star_single_Aligned.sortedByCoord.out.bam" ]; then
    echo "Merging BAM files..."
    samtools merge -@ 16 "$output_dir/final_merged.bam" \
        "$output_dir/star_paired_Aligned.sortedByCoord.out.bam" \
        "$output_dir/star_single_Aligned.sortedByCoord.out.bam"
elif [ -f "$output_dir/star_paired_Aligned.sortedByCoord.out.bam" ]; then
    echo "Using only paired-end BAM for final BAM..."
    cp "$output_dir/star_paired_Aligned.sortedByCoord.out.bam" "$output_dir/final_merged.bam"
elif [ -f "$output_dir/star_single_Aligned.sortedByCoord.out.bam" ]; then
    echo "Using only singleton BAM for final BAM..."
    cp "$output_dir/star_single_Aligned.sortedByCoord.out.bam" "$output_dir/final_merged.bam"
else
    echo "ERROR: No valid BAM files found!"
    exit 1
fi

# 5️⃣ Index final BAM
echo "Indexing final BAM..."
samtools index "$output_dir/final_merged.bam"

echo "✅ STAR Alignment and BAM merging completed successfully!"

