#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sc-blat.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sc-blat.err
### run this script after sc_fusions.ipynb and before snRNA-final.ipynb ###
### define directories ###
BLAT_EXEC="/nfs/users/nfs_l/lr26/blat"
GENOME_FILE="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
FASTA_DIR="/lustre/scratch126/cellgen/behjati/lr26/snRNA/"
OUTPUT_FILE="/lustre/scratch126/cellgen/behjati/lr26/snRNA/blat_top5-parallel.csv"
TEMP_DIR="/lustre/scratch126/cellgen/behjati/lr26/snRNA/temp"

### make the temporary directory ###
mkdir -p "$TEMP_DIR"

### change to home ###
cd /nfs/users/nfs_l/lr26/

### create output file header ###
echo "Sequence,Matches,MisMatches,RepMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts" > "$OUTPUT_FILE"

### run blat on a single sequence creating a psl file ###
run_blat() {
    single_fasta="$1"
    sequence_name=$(basename "$single_fasta" .fasta)
    psl_file="$TEMP_DIR/output_${sequence_name}.psl"
    ### Run BLAT on each of the sequences ###
    "$BLAT_EXEC" "$GENOME_FILE" "$single_fasta" "$psl_file" -out=psl -minIdentity=70
    ### process the psl to extract top 5 hits per sequence ###
    tail -n +6 "$psl_file" | LC_ALL=C sort -k10,10 -k1,1nr | \
    awk -F'\t' '{
        count[$10]++
        if (count[$10] <= 5) print
    }' | while IFS=$'\t' read -r matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts; do
        echo "$sequence_name,$matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts" >> "$OUTPUT_FILE"
    done
}

### then loop over each of the fastas in the directory ###
for fasta_file in "$FASTA_DIR"/*.fasta; do
    file_prefix=$(basename "$fasta_file" .fasta)
    ### split into single-sequence fasta files ###
    csplit -z -f "$TEMP_DIR/${file_prefix}_" -b "%03d.fasta" "$fasta_file" '/^>/' '{*}' > /dev/null
    ### process each sequence individually ###
    for seq_fasta in "$TEMP_DIR/${file_prefix}_"*.fasta; do
        run_blat "$seq_fasta"
        rm -f "$seq_fasta"
    done
done

rm -rf "$TEMP_DIR"
