#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-blast-mom.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-blast-mom.err
set -euo pipefail

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate my-python

cd /lustre/scratch126/cellgen/behjati/lr26/blast-search

#!/bin/bash

# === Configurable paths ===
QUERY="/lustre/scratch126/cellgen/behjati/lr26/blast-search/d4z4.fasta"  # D4Z4 consensus
TARGET_CONTIGS="/nfs/team274/lr26/PacBio-mom/PacBio-mom.bp.p_ctg.fasta"  # Patient's contigs
CONTIG_BAM="/nfs/team274/lr26/PacBio-mom/PacBio-mom.bp.p_ctg_vs_t2t.bam"  # BAM of contigs aligned to T2T

DB_PREFIX="mom_genome_db"
THREADS=8
PIDENTITY_THRESHOLD=80
QCOV_THRESHOLD=85

# === Output files ===
BLAST_OUTPUT="d4z4_matches_mom.tsv"
FILTERED_HITS="d4z4_confident_mom.tsv"
REGIONS_BED="d4z4_regions_mom.txt"
EXTRACTED_FASTA="extracted_d4z4s_mom.fasta"
D4Z4_DB="d4z4_sequences_mom_db"
ALL_VS_ALL_OUT="d4z4_sequences_all_vs_all_mom.tsv"
MAPPED_COORDS="d4z4_mom_regions_mapped.tsv"

# Step 1: Make BLAST DB from patient contigs
echo "Step 1: Create BLAST DB from patient contigs..."
makeblastdb -in "$TARGET_CONTIGS" -dbtype nucl -out "$DB_PREFIX"

# Step 2: Run BLAST to find D4Z4 matches on contigs
echo "Step 2: Run BLAST search for D4Z4 hits on contigs..."
#blastn -query "$QUERY" -db "$DB_PREFIX" \
   # -out "$BLAST_OUTPUT" \
   # -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
   # -num_threads "$THREADS" \
   # -perc_identity "$PIDENTITY_THRESHOLD" \
   # -task blastn \
   # -strand both

# Step 3: Filter BLAST hits for coverage and identity
echo "Step 3: Filter BLAST hits for coverage >= $QCOV_THRESHOLD% and identity >= $PIDENTITY_THRESHOLD%..."
#QLEN=$(grep -v "^>" "$QUERY" | tr -d '\n' | wc -c)
#awk -v qlen="$QLEN" '$3 >= 80 && $4 >= (0.85 * qlen)' "$BLAST_OUTPUT" > "$FILTERED_HITS"


# Step 4: Extract contig coordinates of confident hits
echo "Step 4: Extract contig coordinates from filtered hits..."
#awk '{if ($9 < $10) print $2 ":" $9 "-" $10; else print $2 ":" $10 "-" $9}' "$FILTERED_HITS" > "$REGIONS_BED"

# Step 5: Extract D4Z4 sequences from contigs
echo "Step 5: Extract D4Z4 sequences from contigs..."
#xargs samtools faidx "$TARGET_CONTIGS" < "$REGIONS_BED" > "$EXTRACTED_FASTA"

# *** NEW STEP: Map contig D4Z4 regions to T2T chromosomes using BAM alignments ***
# === Step 6: Map contig D4Z4 regions to T2T chromosomes using BAM alignments (Python) ===
echo "Step 6: Map contig D4Z4 regions to T2T reference coordinates using Python..."

python3 /nfs/users/nfs_l/lr26/shells/map_contigs_to_T2T.py \
    "$CONTIG_BAM" \
    "$REGIONS_BED" \
    "$MAPPED_COORDS"

# Rename FASTA headers in extracted FASTA using mapped coordinates
echo "Renaming FASTA headers with T2T coordinates..."
paste <(grep "^>" "$EXTRACTED_FASTA" | sed 's/>//') <(cut -f2 "$MAPPED_COORDS") | \
while read -r original new; do
    sed -i "s|>$original|>$new|" "$EXTRACTED_FASTA"
done

# Step 7: Build BLAST DB of extracted D4Z4 sequences
echo "Step 7: Build BLAST DB from extracted D4Z4 sequences..."
makeblastdb -in "$EXTRACTED_FASTA" -dbtype nucl -out "$D4Z4_DB"

# Step 8: Run all-vs-all BLAST of extracted D4Z4 sequences
echo "Step 8: Run all-vs-all BLAST search..."
blastn -query "$EXTRACTED_FASTA" -db "$D4Z4_DB" \
  -out "$ALL_VS_ALL_OUT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -perc_identity 30 -dust no -soft_masking false -ungapped -max_hsps 1 -strand both \
  -num_threads "$THREADS"
