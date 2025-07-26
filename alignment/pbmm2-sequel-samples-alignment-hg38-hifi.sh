#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-alignment-hg38-lucie_job
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/hg38.mmi"
input_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio_raw"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-hg38-hifi"

# List of subdirectories to process
sub_dirs=("1_A01" "1_A02" "2_B01")  # Add more as needed

# Process each subdirectory
for sub_dir in "${sub_dirs[@]}"; do
    full_input_dir="$input_dir/$sub_dir"
    full_output_dir="$output_dir/$sub_dir"
    tmp_dir="$full_output_dir/tmp"

    # Check if input directory exists
    if [ ! -d "$full_input_dir" ]; then
        echo "Warning: Input directory $full_input_dir does not exist. Skipping..."
        continue
    fi

    # Create output and temp directories if they don’t exist
    mkdir -p "$full_output_dir" "$tmp_dir"

    # Set TMPDIR for pbmm2 (samtools sorting)
    export TMPDIR="$tmp_dir"

    # Process each BAM file in the subdirectory
    for bam_file in "$full_input_dir"/*.bam; do
        if [ -f "$bam_file" ]; then
            base_name=$(basename "$bam_file" .bam)
            output_bam="$full_output_dir/${base_name}_pbmm2-farm22-bam.bam"

            echo "Aligning $bam_file to $reference..."
            pbmm2 align "$reference" "$bam_file" "$output_bam" --preset HIFI --sort -j 16

            echo "pbmm2 alignment completed: $output_bam"
        else
            echo "Warning: No BAM files found in $full_input_dir."
        fi
    done
done

echo "All alignments completed."
