#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/giabbamfastq.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/giabbamfastq.e
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
input_dir="/lustre/scratch126/casm/team274sb/ExternalData/byUser/na15/PacBio_GIAB_trio/bams/"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/fastq-revio-sprq"

# Define the three specific sample directory names
samples=("m84039_241001_220042_s2.hifi_reads.bc2018.bam" "m84039_241002_000337_s3.hifi_reads.bc2020.bam" "m84039_241002_020632_s4.hifi_reads.bc2021.bam")

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
