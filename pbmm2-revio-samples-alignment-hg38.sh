#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/hg38tum.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/hg38tum.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J pbmm2-alignment-hg38-revio
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.mmi"
input_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio_Revio_raw"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio-nontared"

sub_dirs=("1_B01")

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

    # Process only BAM files that contain "hifi_reads" in their name
    for bam_file in "$full_input_dir"/*hifi_reads*.bam; do
        if [ -f "$bam_file" ]; then
            base_name=$(basename "$bam_file" .bam)
            output_bam="$full_output_dir/${base_name}_pbmm2-farm22-bam.bam"

            echo "Aligning $bam_file to $reference..."
            pbmm2 align "$reference" "$bam_file" "$output_bam" --preset CCS --sort -j 16 --unmapped

            echo "pbmm2 alignment completed: $output_bam"
        fi
    done
done

echo "All alignments completed."
