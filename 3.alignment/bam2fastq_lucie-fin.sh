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

### Directories to process
input_dir="/lustre/scratch126/casm/team274sb/na15/PacBio/raw_output"
output_dir="/lustre/scratch126/casm/team274sb/lr26/fastq"

for sub_dir in "1_A01" "2_B01"; do
  full_dir="$input_dir/$sub_dir"
  if [ -d "$full_dir" ]; then
    for bam_file in "$full_dir"/*.bam; do
      if [ -f "$bam_file" ]; then
        base_name=$(basename "$bam_file" .bam)
        bam2fastq "$bam_file" -o "$output_dir/$base_name.fastq"
        gzip "$output_dir/$base_name.fastq"
      fi
    done
  else
    echo "Directory $full_dir does not exist. Skipping."
  fi
done
