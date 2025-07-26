#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q normal
#BSUB -J ngmlr-alignment-t2t-lucie_job
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference_gz="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa.gz"
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/na15/PacBio/raw_output"
output_dir="/lustre/scratch126/casm/team274sb/lr26/ngmlr-alignment"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Decompress reference if not already done
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi

# Process BAM files
for sub_dir in "1_A01" "2_B01"; do
  full_dir="$input_dir/$sub_dir"

  if [ -d "$full_dir" ]; then
    for bam_file in "$full_dir"/*.bam; do
      if [ -f "$bam_file" ]; then
        base_name=$(basename "$bam_file" .bam)
        output_bam="$output_dir/${base_name}_ngmlr.bam"

        echo "Aligning $bam_file to $reference using NGMLR..."
        ngmlr -r "$reference" -q "$bam_file" -o "$output_bam" -t 16 -x pacbio

      fi
    done
  else
    echo "Warning: Directory $full_dir does not exist. Skipping."
  fi
done

echo "NGMLR alignment completed."
