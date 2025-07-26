#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 10
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J bam2fastq_job
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories
input_dir="/lustre/scratch126/casm/team274sb/lr26/Revio_raw/"
output_dir="/lustre/scratch126/casm/team274sb/lr26/fastq-revio"

### Make sure the output directory exists
mkdir -p "$output_dir"

for sub_dir in "$input_dir"/*; do  # Remove the trailing slash in the glob
  sub_dir=$(basename "$sub_dir")  # Extract directory name only
  full_dir="$input_dir/$sub_dir"

  if [ -d "$full_dir" ]; then
    echo "Processing directory: $full_dir"

    bam_files=("$full_dir"/*.bam)
    if [ ! -e "${bam_files[0]}" ]; then
      echo "No BAM files found in $full_dir. Skipping."
      continue
    fi

    for bam_file in "${bam_files[@]}"; do
      base_name=$(basename "$bam_file" .bam)
      
      # Convert BAM to FASTQ
      bam2fastq "$bam_file" -o "$output_dir/$base_name"

      # Check if the output file was created before compressing
      if [ -f "$output_dir/$base_name.fastq" ]; then
        if command -v pigz &> /dev/null; then
          pigz "$output_dir/$base_name.fastq"
        else
          gzip "$output_dir/$base_name.fastq"
        fi
      else
        echo "Warning: FASTQ file for $bam_file was not created."
      fi
    done
  fi
done
