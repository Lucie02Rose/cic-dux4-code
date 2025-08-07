#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-blast-t2t.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-blast-t2t.err

### failure management ###
set -euo pipefail

### paths and variables ###
### consensus, reference, database ###
d4z4="/lustre/scratch126/casm/team274sb/lr26/T2T/d4z4.fasta"  ### found here in the 6.d4z4_structural_protein folder           
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"                         
database="target_genome_db"    

### outputs ###
blast_output="/lustre/scratch126/casm/team274sb/lr26/T2T/d4z4_matches.tsv"
filtered_d4z4s="/lustre/scratch126/casm/team274sb/lr26/T2T/d4z4_confident.tsv"
d4z4_regions="/lustre/scratch126/casm/team274sb/lr26/T2T/d4z4_regions.txt"
extracted_fastas="/lustre/scratch126/casm/team274sb/lr26/T2T/extracted_d4z4s.fasta"
d4z4_db="d4z4_sequences_db"
output_tsv="/lustre/scratch126/casm/team274sb/lr26/T2T/d4z4_sequences_all_vs_all.tsv"

### change to the correct directory ###
cd /lustre/scratch126/casm/team274sb/lr26/T2T/

### make blast database from genome ###
echo "make db from genome"
makeblastdb -in "$reference" -dbtype nucl -out "$database"

### make blast database from genome ###
echo "run blast of d4z4 against the genome"
blastn -query "$d4z4" -db "$database" \
    -out "$blast_output" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads 8 \
    -perc_identity 85 \
    -task blastn \
    -strand both

### get query length and filter based on coverage and identity ###
echo "filter blasts for full coverage (95%) and identity (85%)"
QLEN=$(grep -v "^>" "$d4z4" | tr -d '\n' | wc -c)
awk -v qlen="$QLEN" '$3 >= 85 && $4 >= (0.95 * qlen)' "$blast_output" > "$filtered_d4z4s"

### extract the coodinates for confident sequences ###
echo "extract coordinates for confident hits"
awk '{if ($9 < $10) print $2 ":" $9 "-" $10; else print $2 ":" $10 "-" $9}' "$filtered_d4z4s" > "$d4z4_regions"

### extract those regions into corresponding fastas ###
echo "extract D4Z4 sequences from genome"
xargs samtools faidx "$reference" < "$d4z4_regions" > "$extracted_fastas"

### build a blast database of the extracted d4z4s ###
echo "blast db of extracted sequences"
makeblastdb -in "$extracted_fastas" -dbtype nucl -out "$d4z4_db"

### run an all-by-all blast making sure that even the low hits kept ###
echo "all-vs-all blast"
blastn -query "$extracted_fastas" -db "$d4z4_db" \
    -out "$output_tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 10 \
    -num_threads 8 \
    -dust no \
    -soft_masking false \
    -ungapped \
    -max_hsps 1 \
    -perc_identity 30 \
    -task blastn \
    -strand both

### print statement - completed ###
echo "completed"
