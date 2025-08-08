#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q normal
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J.err

### activate the base conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### fasta files of soft clipped reads ###
fasta_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-aligned/fasta"
### blastn database made from the T2T reference fasta file ###
blast_db="/lustre/scratch126/casm/team274sb/lr26/T2T/T2T_chm13_db"
### output directory for blast results ###
output_dir="/lustre/scratch126/casm/team274sb/lr26/PacBio-aligned/blast_results"
#### csv file to store all blast results ###
csv_output="/lustre/scratch126/casm/team274sb/lr26/PacBio-aligned/blast_results/all_blast_results.csv"

### create the output directory ###
mkdir -p "$output_dir"

### make the csv file header ###
echo "Read_ID,Chromosome,Percent_Identity,Alignment_Length,Mismatches,Gap_Openings,Start,End,Query_Start,Query_End,E_value,Bit_Score" > "$csv_output"

### loop through each read in the fasta directory ###
for fasta_file in "$fasta_dir"/*.fasta; do
    ### base name of each sequence id ###
    base_name=$(basename "$fasta_file" .fasta)
    ### define output per sequence id ###
    output_file="$output_dir/${base_name}_blast_results.txt"
    ### run blastn and save the results ###
    raw_output=$(blastn -query "$fasta_file" -db "$blast_db" -outfmt "6" -max_target_seqs 5)
    ### sort output by bit score - best to worse, retaining top 5 ###
    top_5_results=$(echo "$raw_output" | sort -k12,12gr | head -n 5)
    ### append to the final csv file ###
    ### loop through each and append the top 5 ###
    while IFS=$'\t' read -r read_id chrom percent_identity alignment_length mismatches gap_openings start end query_start query_end e_value bit_score; do
        echo "$read_id,$chrom,$percent_identity,$alignment_length,$mismatches,$gap_openings,$start,$end,$query_start,$query_end,$e_value,$bit_score" >> "$csv_output"
    done <<< "$top_5_results"
    ### check if the output is empty and print if it is ###
    if [ ! -s "$output_file" ]; then
        echo "No BLAST results for $fasta_file"
    fi
done
### print statement for being finished ###
echo "BLAST search completed for all FASTA files."
