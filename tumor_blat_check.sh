#!/bin/bash
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Input file
fasta_file="/lustre/scratch126/casm/tream274sb/lr26/filtered_bams/tumor.fasta"
output_file="/lustre/scratch126/casm/tream274sb/lr26/filtered_bams/tumor_regions_interest.fasta"

# Define the range of interest
range_start=128410000
range_end=128423000

# Create or clear the output file
> "$output_file"

# Loop through each sequence in the fasta file
while read -r header; do
    # Extract the start and end coordinates from the header
    if [[ $header =~ start([0-9]+).*end([0-9]+) ]]; then
        start=${BASH_REMATCH[1]}
        end=${BASH_REMATCH[2]}

        # Check if the start and end overlap with the range of interest
        if [[ ($start -ge $range_start && $start -le $range_end) || 
              ($end -ge $range_start && $end -le $range_end) || 
              ($start -le $range_start && $end -ge $range_end) ]]; then

            # Append the header to the output file
            echo "$header" >> "$output_file"
            
            # Read and append the sequence (the next line after header)
            read sequence
            echo "$sequence" >> "$output_file"
        fi
    else
        # Append the header if it doesn't contain start/end coordinates
        echo "$header" >> "$output_file"
        read sequence
        echo "$sequence" >> "$output_file"
    fi
done < "$fasta_file"

echo "Sequences extracted and saved to $output_file"

