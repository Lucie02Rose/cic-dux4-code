#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Define paths
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"   # Update if BLAT is elsewhere
GENOME_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"  # Genome reference
FASTA_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/fasta/"  # FASTA input files directory
OUTPUT_FILE="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/blat_top5.csv"  # Final output file
TEMP_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/temp"  # Use RAM for temporary processing

# Ensure temp directory exists
mkdir -p "$TEMP_DIR"

# Create/clear the output file and add CSV header
echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$OUTPUT_FILE"

# Function to run BLAT for a single FASTA file
run_blat() {
    fasta_file="$1"
    sequence_name=$(basename "$fasta_file" .fasta)  # Extract file name

    # Run BLAT and store PSL output in RAM
    "$BLAT_EXEC" "$GENOME_FILE" "$fasta_file" "$TEMP_DIR/output_${sequence_name}.psl" -out=ps -minIdentity=90

    # Process PSL output: Remove header, sort by matches, get top 5, store in memory
    results=""
    while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
        results+="$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts"$'\n'
    done < <(tail -n +6 "$TEMP_DIR/output_${sequence_name}.psl" | LC_ALL=C sort -k1,1nr | head -5)

    # Append results to final CSV file
    echo -n "$results" >> "$OUTPUT_FILE"

    # Cleanup
    rm -f "$TEMP_DIR/output_${sequence_name}.psl"
}

# Process each FASTA file sequentially
for fasta_file in "$FASTA_DIR"/*.fasta; do
    run_blat "$fasta_file"
done

# Cleanup temp directory
rm -rf "$TEMP_DIR"
