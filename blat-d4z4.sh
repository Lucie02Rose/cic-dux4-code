#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

# Define paths
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"   # Path to BLAT executable
GENOME_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"  # Genome reference
FASTA="/lustre/scratch126/casm/team274sb/lr26/T2T/d4d4_repeat_sequence-chr4.fasta"  # Input FASTA file
OUTPUT_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_d4z4.csv"  # Final output file
TEMP_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/temp_output.psl"  # Temporary PSL output

# Create/clear the output file and add CSV header
echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$OUTPUT_FILE"

# Run BLAT with minIdentity=85 and output to PSL file
"$BLAT_EXEC" -minIdentity=85 "$GENOME_FILE" "$FASTA" "$TEMP_FILE" -out=psl

# Extract the sequence name from the FASTA file name
sequence_name=$(basename "$FASTA" .fasta)

# Process PSL output: Remove header, sort by matches, and append to CSV
grep -v '^#' "$TEMP_FILE" | LC_ALL=C sort --parallel=16 -k1,1nr | while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
    echo "$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts"
done >> "$OUTPUT_FILE"

# Cleanup temporary PSL file
rm -f "$TEMP_FILE"

echo "BLAT search completed. Results saved in $OUTPUT_FILE"
