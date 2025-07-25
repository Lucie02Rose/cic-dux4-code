#BSUB -o /nfs/users/nfs_l/lr26/outputs/giabbamfastq-%J.o
#BSUB -e /nfs/users/nfs_l/lr26/errors/giabbamfastq-%J.e
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J bam2fastq_job
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq"

# Define the three specific sample directory names
samples=("blood_1C01_hifi_reads.bam" "tumor_1A01_hifi_reads.bam" "tumor_1A02_hifi_reads.bam" "tumor_1B01_hifi_reads.bam" "tumor_2B01_hifi_reads.bam" "mom_1B02_hifi_reads.bam")

mkdir -p "$output_dir"

for bam_file_name in "${samples[@]}"; do
  bam_file="$input_dir/$bam_file_name"

  if [ -f "$bam_file" ]; then
    echo "Processing file: $bam_file"
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
  else
    echo "File $bam_file does not exist. Skipping."
  fi
done
