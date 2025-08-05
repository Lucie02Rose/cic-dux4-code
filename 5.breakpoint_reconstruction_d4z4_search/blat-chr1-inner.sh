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
FASTA_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/fasta-chr1-inner/"  # FASTA input files directory
OUTPUT_FILE="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/blat_top5-chr1-inner-parallel.csv"  # Final output file
TEMP_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/temp-chr1"  # Use RAM for temporary processing

# Ensure temp directory exists
mkdir -p "$TEMP_DIR"

# Create/clear the output file and add CSV header
echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$OUTPUT_FILE"

# Function to run BLAT for a single FASTA file
run_blat() {
    fasta_file="$1"
    sequence_name=$(basename "$fasta_file" .fasta)  # Extract file name

    # Run BLAT and store PSL output in RAM
    "$BLAT_EXEC" "$GENOME_FILE" "$fasta_file" "$TEMP_DIR/output_${sequence_name}.psl" -out=psl -minIdentity=90

    # Process PSL output: Remove header, sort by matches, get top 5, store in memory
    results=""
    while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
        results+="$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts"$'\n'
    done < <(tail -n +6 "$TEMP_DIR/output_${sequence_name}.psl" | LC_ALL=C sort --parallel=16 -k1,1nr | head -5)

    # Append results to final CSV file
    echo -n "$results" >> "$OUTPUT_FILE"

    # Cleanup
    rm -f "$TEMP_DIR/output_${sequence_name}.psl"
}

export -f run_blat
export BLAT_EXEC GENOME_FILE TEMP_DIR OUTPUT_FILE  # Export variables for parallel

# Find all FASTA files and process them in parallel (max 16 jobs at a time)
find "$FASTA_DIR" -name "*.fasta" | parallel -j 16 run_blat

# Cleanup temp directory
rm -rf "$TEMP_DIR"

#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### this is the core script for the BLAT seach of the breakpoints of interest around the CIC region ###
### prior to this script, I have used functions in the pysam_sequence_blatting.ipynb to extract the soft clips from all breakpoints and turn them into fastas ###
### therefore, define the bam file and region of interest (here it is start, end = 45055000, 45160000 in the Python script)
### then use the extract_soft_clipped_sequences and save_to_csv, then sanitize_read_id, save_to_fasta and convert_csv_to_fasta.
### the Python script makes a fasta file for each sequence based on read id ###
### thir script then runs a parallel blat ###
### do not really need the conda environment ###
### executable blat programme is found here ###
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"   
### genome reference ###
GENOME_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"  
### fasta files of soft clipped reads ###
FASTA_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/fasta/" 
### top 5 blat hits per soft clipped read ###
OUTPUT_FILE="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/blat_top5-parallel.csv" 
### set up a temporary directory for sorting ###
TEMP_DIR="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-tumor-new/temp"  

### make the temporary directory ###
mkdir -p "$TEMP_DIR"
### make the ouput file header depending on the blat programme specifications ###
echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$OUTPUT_FILE"
### function that runs the blat search on the single fasta file - sequence ###
run_blat() {
    fasta_file="$1"
    sequence_name=$(basename "$fasta_file" .fasta)
    psl_file="$TEMP_DIR/output_${sequence_name}.psl"
    csv_file="$TEMP_DIR/${sequence_name}.csv"
    ### the actual command that runs blat ###
    "$BLAT_EXEC" "$GENOME_FILE" "$fasta_file" "$psl_file" -out=psl -minIdentity=85
    ### this is the header for the per-sequence blat ouput ###
    echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$csv_file"
    ### the psl file created in the blat ouput is sorted according to the top best matches and parsed so that only the top 5 would be retained ###
    ### the top 5 psl are then written into the pre-prepared csv file which already has the header prepared ###
    tail -n +6 "$psl_file" | LC_ALL=C sort -k1,1nr | head -5 | while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
        echo "$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts" >> "$csv_file"
    done
    ### the huge psl is deleted not to take up space as it is no needed anymore ####
    rm -f "$psl_file"
}
### export the function for parallel to recognise ###
export -f run_blat
### export the variables that were defined above for parallel to recognise for each parallelised job ###
export BLAT_EXEC GENOME_FILE TEMP_DIR OUTPUT_FILE
### let parallel seach for the fastas and let them run, 16 at a time each using a cpu ###
find "$FASTA_DIR" -name "*.fasta" | parallel -j 16 run_blat
### cleanup of the temporary directory ###
rm -rf "$TEMP_DIR"
