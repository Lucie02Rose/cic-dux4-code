#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-blast-t2t.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-blast-t2t.err
set -euo pipefail

cd /lustre/scratch126/cellgen/behjati/lr26/T2T/

# === Configurable paths ===
QUERY="/lustre/scratch126/cellgen/behjati/lr26/T2T/d4z4.fasta"                  # D4Z4 consensus sequence
TARGET_GENOME="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"                             # Reference or assembly
DB_PREFIX="target_genome_db"                             # Temporary BLAST DB name
THREADS=8
PIDENTITY_THRESHOLD=85
QCOV_THRESHOLD=95

# === Output files ===
BLAST_OUTPUT="/lustre/scratch126/cellgen/behjati/lr26/T2T/d4z4_matches.tsv"
FILTERED_HITS="/lustre/scratch126/cellgen/behjati/lr26/T2T/d4z4_confident.tsv"
REGIONS_BED="/lustre/scratch126/cellgen/behjati/lr26/T2T/d4z4_regions.txt"
EXTRACTED_FASTA="/lustre/scratch126/cellgen/behjati/lr26/T2T/extracted_d4z4s.fasta"
D4Z4_DB="d4z4_sequences_db"
ALL_VS_ALL_OUT="/lustre/scratch126/cellgen/behjati/lr26/T2T/d4z4_sequences_all_vs_all.tsv"

echo "Step 1: Create BLAST DB from genome..."
makeblastdb -in "$TARGET_GENOME" -dbtype nucl -out "$DB_PREFIX"

echo "Step 2: Run BLAST to find D4Z4 matches..."
blastn -query "$QUERY" -db "$DB_PREFIX" \
    -out "$BLAST_OUTPUT" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads "$THREADS" \
    -perc_identity "$PIDENTITY_THRESHOLD" \
    -task blastn \
    -strand both

echo "Step 3: Filter BLAST hits for full coverage (>=95%) and identity (>=85%)..."
# Get query length
QLEN=$(grep -v "^>" "$QUERY" | tr -d '\n' | wc -c)
awk -v qlen="$QLEN" '$3 >= 85 && $4 >= (0.95 * qlen)' "$BLAST_OUTPUT" > "$FILTERED_HITS"

echo "Step 4: Extract coordinates for confident hits..."
awk '{if ($9 < $10) print $2 ":" $9 "-" $10; else print $2 ":" $10 "-" $9}' "$FILTERED_HITS" > "$REGIONS_BED"

echo "Step 5: Extract D4Z4 sequences from genome..."
xargs samtools faidx "$TARGET_GENOME" < "$REGIONS_BED" > "$EXTRACTED_FASTA"

echo "Step 6: Build BLAST DB of extracted D4Z4 sequences..."
makeblastdb -in "$EXTRACTED_FASTA" -dbtype nucl -out "$D4Z4_DB"

echo "Step 7: Run all-vs-all BLAST..."
blastn -query "$EXTRACTED_FASTA" -db "$D4Z4_DB" \
    -out "$ALL_VS_ALL_OUT" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 10 \
    -num_threads "$THREADS" \
    -dust no \
    -soft_masking false \
    -ungapped \
    -max_hsps 1 \
    -perc_identity 30 \
    -task blastn \
    -strand both

echo "✅ Pipeline complete. Output:"
echo "- Filtered hits: $FILTERED_HITS"
echo "- Extracted D4Z4s: $EXTRACTED_FASTA"
echo "- All-vs-all BLAST: $ALL_VS_ALL_OUT"
