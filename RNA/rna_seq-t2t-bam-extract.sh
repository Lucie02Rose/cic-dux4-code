#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
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
reference="/nfs/users/nfs_l/lr26/PacBio/T2T/chm13v2.0.fa"
output_dir="/lustre/scratch126/casm/team274sb/lr26/rna-t2t"
reference_index="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t/T2T_index"

# Define an array of GSM sample IDs
samples=("GSM5024895" "GSM5024896" "GSM5024897" "GSM5024898" "GSM5024899" "GSM5024900" "GSM5024901" "GSM5024902")

# Define base directory where BAM files are stored
base_dir="/lustre/scratch127/cellgen/cellgeni/tickets/tic-3785/results/STAR-Fusion"

# Create the output directory if it doesn't exist
cd "$output_dir"

process_sample() {
    sample=$1
    bam_file="${base_dir}/${sample}/STAR-Fusion_outdir/Aligned.out.bam"

    r1_out="${output_dir}/${sample}_R1.fastq"
    r2_out="${output_dir}/${sample}_R2.fastq"
    single_out="${output_dir}/${sample}_singletons.fastq"

    # Check if BAM file exists before processing
    if [[ -f "$bam_file" ]]; then
        echo "Processing $sample..."
        samtools fastq -1 "$r1_out" -2 "$r2_out" -s "$single_out" "$bam_file"
        echo "Finished extracting FASTQ for $sample"
    else
        echo "Warning: BAM file for $sample not found, skipping..."
    fi
}

# Process samples with controlled concurrency (max 16 jobs at a time)
for sample in "${samples[@]}"; do
    process_sample "$sample" &
    
    # Limit concurrent jobs to 16
    if (( $(jobs | wc -l) >= 16 )); then
        wait -n  # Wait for at least one job to finish before starting a new one
    fi
done
wait  # Ensure all background jobs finish before moving on

