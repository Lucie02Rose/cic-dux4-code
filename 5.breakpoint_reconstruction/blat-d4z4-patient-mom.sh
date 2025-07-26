#!/bin/bash
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

# Define paths
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"
GENOME_FILE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
FASTA_DIR="/lustre/scratch126/casm/team274sb/lr26/filtered_bams"
D4Z4_REPEAT="/lustre/scratch126/casm/team274sb/lr26/T2T/d4d4_repeat_sequence-chr4.fasta"
OUTPUT_DIR="/lustre/scratch126/casm/team274sb/lr26/filtered_bams/blat_mom_patient_results/"
MOM="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/mom.psl"
TUMOR="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/tumor.psl"
BLOOD="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/blood.psl"
TSV_MOM="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/mom.tsv"
TSV_TUMOR="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/tumor.tsv"
TSV_BLOOD="/lustre/scratch126/casm/team274sb/lr26/T2T/blat_results/blood.tsv"

# BAM files
#declare -A BAM_FILES=(
 #   ["tumor"]="tumor-all-chr1410.bam"
  #  ["blood"]="patient-blood-chr1410.bam"
   # ["mom"]="mom-chr1410.bam"
#)

# Ensure output directory exists
#mkdir -p "$OUTPUT_DIR"

# Step 1: Convert BAM to FASTA using bedtools + seqtk
#for sample in "${!BAM_FILES[@]}"; do
 #   BAM_PATH="$FASTA_DIR/${BAM_FILES[$sample]}"
  #  FASTA_OUT="$OUTPUT_DIR/${sample}.fasta"

   # echo "Extracting FASTA from $BAM_PATH..."
    #samtools fasta "$BAM_PATH" > "$FASTA_OUT"

#done


# Function to process a single PSL file
process_psl() {
    PSL_FILE=$1
    TSV_FILE=$2
    SAMPLE_NAME=$3  # Pass sample name for logging

    echo "Processing BLAT results for $SAMPLE_NAME..."
    echo -e "Sample\tMatches\tMisMatches\tRepMatches\tnCount\tqNumInsert\tqBaseInsert\ttNumInsert\ttBaseInsert\tstrand\tqName\tqSize\tqStart\tqEnd\ttName\ttSize\ttStart\ttEnd\tblockCount\tblockSizes\tqStarts\ttStarts" > "$TSV_FILE"

    awk 'NR > 5 {print "'"$SAMPLE_NAME"'\t"$0}' "$PSL_FILE" | sort --parallel=16 -k1,1nr >> "$TSV_FILE"

    # Cleanup PSL file (only if successful)
    if [ $? -eq 0 ]; then
        rm -f "$PSL_FILE"
    else
        echo "Error processing $SAMPLE_NAME, PSL file retained for debugging."
    fi
}

# Run the function for each file
process_psl "$MOM" "$TSV_MOM" "MOM"
process_psl "$TUMOR" "$TSV_TUMOR" "TUMOR"
process_psl "$BLOOD" "$TSV_BLOOD" "BLOOD"
