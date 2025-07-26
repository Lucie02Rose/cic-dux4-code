#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J ngmlr-sort-index-t2t
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

#!/bin/bash

# Define the base directory
BASE_DIR="/lustre/scratch126/casm/team274sb/lr26/ngmlr-alignment"

# Check if the directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Error: Directory $BASE_DIR does not exist."
    exit 1
fi

# Loop through each subdirectory
for subdir in "$BASE_DIR"/*/; do
    # Ensure it's a directory
    if [ -d "$subdir" ]; then
        echo "Processing directory: $subdir"
        
        # Find all BAM files in this subdirectory
        find "$subdir" -maxdepth 1 -type f -name "*.bam" | while read -r bam; do
            echo "Processing BAM file: $bam"
            
            # Define the sorted BAM file name
            sorted_bam="${bam%.bam}.sorted.bam"
            
            # Sort BAM file
            samtools sort -o "$sorted_bam" "$bam"
            
            # Check if sorting was successful
            if [ $? -eq 0 ]; then
                echo "Sorted: $sorted_bam"
                
                # Index the sorted BAM file
                samtools index "$sorted_bam"
                
                if [ $? -eq 0 ]; then
                    echo "Indexed: ${sorted_bam}.bai"
                else
                    echo "Error: Failed to index $sorted_bam"
                fi
            else
                echo "Error: Failed to sort $bam"
            fi
        done
    fi
done

echo "All BAM files processed."
