#!/bin/bash
#BSUB -n 64
#BSUB -M 350000
#BSUB -R 'span[hosts=1] select[mem>350000] rusage[mem=350000]'
#BSUB -q hugemem
#BSUB -J blat_top20
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Define paths
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"  # Update if BLAT is elsewhere
GENOME_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"  # Reference genome
TEMP_DIR="/lustre/scratch126/casm/team274sb/lr26/rolling_window/temp"  # Folder with all seq_x.fasta files
OUTPUT_FILE="/lustre/scratch126/casm/team274sb/lr26/rolling_window/blat_top20.csv"  # Final output

# Ensure temp directory exists
cd "$TEMP_DIR"

# Clear previous output and add CSV header
echo "Query,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,BlockCount,BlockSizes,qStarts,tStarts" > "$OUTPUT_FILE"

# Function to run BLAT for a single sequence file
run_blat() {
    fasta_file="$1"
    sequence_name=$(basename "$fasta_file" .fasta)  # Extract sequence name

    echo "Running BLAT for $sequence_name..."

    # Run BLAT and save PSL output
    "$BLAT_EXEC" "$GENOME_FILE" "$fasta_file" "$TEMP_DIR/output_${sequence_name}.psl" -out=psl -tileSize=7 -minMatch=1 -oneOff=1 -repMatch=1000000 -maxIntron=1000000 -stepSize=3 -minScore=10

    # Optional: Add a header if needed, can be skipped since BLAT already adds it
    echo "# BLAT output for $sequence_name" > "$TEMP_DIR/output_${sequence_name}_with_header.psl"
    cat "$TEMP_DIR/output_${sequence_name}.psl" >> "$TEMP_DIR/output_${sequence_name}_with_header.psl"

    # Extract top 20 hits and append to CSV
    results=""
    while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
        results+="$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts\n"
    done < <(tail -n +6 "$TEMP_DIR/output_${sequence_name}_with_header.psl" | LC_ALL=C sort --parallel=16 -k1,1nr | head -20)

    # Append results only if there is data
    if [[ -n "$results" ]]; then
        echo -e "$results" >> "$OUTPUT_FILE"
    fi

    # Optional: Keep the PSL file with header for further reference
    # rm -f "$TEMP_DIR/output_${sequence_name}.psl"
}


export -f run_blat
export BLAT_EXEC GENOME_FILE TEMP_DIR OUTPUT_FILE  # Export variables for parallel

# Step 1: Process all individual sequences in parallel
find "$TEMP_DIR" -name "*.fasta" | parallel -j 16 run_blat

echo "BLAT processing complete. Results saved in: $OUTPUT_FILE"
