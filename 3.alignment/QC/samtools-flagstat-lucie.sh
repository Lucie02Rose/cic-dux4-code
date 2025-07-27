#!/bin/bash
### define parameters for the LSF ###
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q long
#BSUB -J flagstat
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### activate conda environment - samtools is in base ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define the directories ###
pbmm2_output="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment"

### this searches for all aligned bam diles in the subdirectories on pbmm2-alignment
### it loops through and checks if it is a directory or if there is a valid bam file
### if it finds a bamfile then it defines a flagstat report incorporating the bam file name
### runs a multithreaded samtools flagstat report which outputs some alignment metrics
for sub_dir in "$pbmm2_output"/*; do
    if [ -d "$sub_dir" ]; then 
        for bam_file in "$sub_dir"/*.bam; do
            if [ -f "$bam_file" ]; then  
              
                flagstat_report="${bam_file%.bam}.flagstat"
                
 
                echo "Running flagstat for $bam_file..."
                samtools flagstat -@ 16 "$bam_file" > "$flagstat_report"
                
                echo "Flagstat report saved: $flagstat_report"
            fi
        done
    fi
done

echo "All flagstat reports generated."
