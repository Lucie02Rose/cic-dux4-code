#!/bin/bash
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q long
#BSUB -J flagstat
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories

pbmm2_output="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38"

# Loop through all subdirectories
for sub_dir in "$pbmm2_output"/*; do
    if [ -d "$sub_dir" ]; then  # Check if it's a directory
        for bam_file in "$sub_dir"/*.bam; do
            if [ -f "$bam_file" ]; then  # Check if BAM file exists
                
                # Extract the filename without the path
                bam_basename=$(basename "$bam_file" .bam)

                echo "Processing: $bam_file"

                # Define output filenames
                filtered_bam="$sub_dir/${bam_basename}.filtered.bam"
                fixed_bam="$sub_dir/${bam_basename}.fixed.bam"
                softclipped_bam="$sub_dir/${bam_basename}.softclipped.bam"

                # Remove supplementary and secondary alignments
                samtools view -F 2308 -b -o "$filtered_bam" "$bam_file"

                # Fix flags and mate information
                samtools fixmate -r "$filtered_bam" "$fixed_bam"

                # Soft-clip long insertions and deletions
                bamtools convert -format sam -in "$fixed_bam" | sed 's/\([0-9]\{3,\}\)[ID]/\1S/g' | samtools view -b -o "$softclipped_bam"

                # Index the final BAM file
                samtools index "$softclipped_bam"

                echo "Finished processing: $softclipped_bam"
            fi
        done
    fi
done

echo "All BAM files processed."
