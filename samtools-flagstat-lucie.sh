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

pbmm2_output="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment"

# Loop through all subdirectories
for sub_dir in "$pbmm2_output"/*; do
    if [ -d "$sub_dir" ]; then  # Check if it's a directory
        for bam_file in "$sub_dir"/*.bam; do
            if [ -f "$bam_file" ]; then  # Check if it's a BAM file
                # Define output flagstat report filename
                flagstat_report="${bam_file%.bam}.flagstat"
                
                # Run samtools flagstat
                echo "Running flagstat for $bam_file..."
                samtools flagstat -@ 16 "$bam_file" > "$flagstat_report"
                
                echo "Flagstat report saved: $flagstat_report"
            fi
        done
    fi
done

echo "All flagstat reports generated."
