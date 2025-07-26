#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q long
#BSUB -J bam-indexing
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Define directories
input_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"
output_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"

# Loop through all BAM files and index them
for bam_file in "$input_dir"/*_Aligned.sortedByCoord.bam; do
    # Ensure the BAM file exists
    if [[ -f "$bam_file" ]]; then
        echo "Indexing BAM file: $bam_file"
        samtools index "$bam_file"
    else
        echo "WARNING: BAM file $bam_file not found!"
    fi
done

echo "BAM indexing completed!"
