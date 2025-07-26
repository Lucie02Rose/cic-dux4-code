#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 300000
#BSUB -R 'span[hosts=1] select[mem>300000] rusage[mem=300000]'
#BSUB -q hugemem
#BSUB -J rna-t2t-alignemnt
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### Directories to process
output_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"

# Move to the output directory
cd "$output_dir"

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

