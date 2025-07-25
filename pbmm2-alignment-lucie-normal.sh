#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q normal
#BSUB -J pbmm2-normal-alignment-t2t-lucie_job
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories to process
reference_gz="/nfs/users/nfs_l/lr26/PacBio/T2T/chm13v2.0.fa.gz"
reference="/nfs/users/nfs_l/lr26/PacBio/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/na15/PacBio/raw_output"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment"

# Decompress reference if not already done
if [ ! -f "$reference" ]; then
  gunzip -c "$reference_gz" > "$reference"
fi

# Indexing the T2T reference genome (if index files don't exist)
if [ ! -f "$reference.mmi" ]; then
  pbmm2 index "$reference"
fi

# Process BAM files
for sub_dir in "1_A01" "2_B01"; do
  full_dir="$input_dir/$sub_dir"
  
  if [ -d "$full_dir" ]; then
    for bam_file in "$full_dir"/*.bam; do
      if [ -f "$bam_file" ]; then
        base_name=$(basename "$bam_file" .bam)
        output_bam="$output_dir/${base_name}_aligned.bam"

        echo "Aligning $bam_file to $reference..."
        pbmm2 align "$reference" "$bam_file" "$output_bam" --preset CCS --sort
      fi
    done
  else
    echo "Warning: Directory $full_dir does not exist. Skipping."
  fi
done

echo "Alignment completed."
