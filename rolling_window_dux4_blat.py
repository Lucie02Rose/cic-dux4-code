#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#rolling window blat
import pysam
import subprocess

# Extract 20bp rolling windows from input sequence
def extract_rolling_windows(input_fasta, window_size=20):
    samfile = pysam.FastaFile(input_fasta)
    windows = []
    for seq_name in samfile.references:
        sequence = samfile.fetch(seq_name)
        for i in range(len(sequence) - window_size + 1):
            windows.append(sequence[i:i+window_size])
    samfile.close()
    return windows

# Save rolling windows to a FASTA file
def save_to_fasta(windows, output_fasta):
    with open(output_fasta, 'w') as f:
        for i, seq in enumerate(windows):
            f.write(f">seq_{i+1}\n{seq}\n")

# Run BLAT for each 20bp sequence and capture output
def run_blat(query_fasta, reference_db, output_psl):
    command = ['/nfs/users/nfs_l/lr26/blat', reference_db, query_fasta, output_psl]
    subprocess.run(command, check=True)

# Parse the BLAT PSL output and get the top 20 hits for each query
def parse_blat_psl(psl_file, top_n=20):
    hits = {}
    
    with open(psl_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.startswith("target"):
                continue  # Skip header lines
            fields = line.strip().split("\t")
            query_name = fields[9]
            score = int(fields[0])  # BLAT score (higher means better)
            
            if query_name not in hits:
                hits[query_name] = []
            
            # Add the hit to the list of hits for this query
            hits[query_name].append((score, line.strip()))
    
    # Sort each query's hits by score and take the top N
    top_hits = {}
    for query_name, hit_list in hits.items():
        # Sort hits by score (descending order)
        hit_list.sort(reverse=True, key=lambda x: x[0])
        top_hits[query_name] = hit_list[:top_n]
    
    return top_hits

# Save top 20 hits for each query to a new file
def save_top_hits(top_hits, output_file):
    with open(output_file, 'w') as f:
        for query_name, hits in top_hits.items():
            f.write(f">Hits_for_{query_name}\n")
            for _, hit in hits:
                f.write(f"{hit}\n")

# Full process combining the steps
def process_fasta(input_fasta, reference_db, output_fasta, output_psl, output_hits_file, window_size=20, top_n=20):
    # Step 1: Extract rolling windows
    windows = extract_rolling_windows(input_fasta, window_size)
    
    # Step 2: Save windows to FASTA
    save_to_fasta(windows, output_fasta)
    
    # Step 3: Run BLAT against the reference database
    run_blat(output_fasta, reference_db, output_psl)
    
    # Step 4: Parse the BLAT output to get the top N hits for each query
    top_hits = parse_blat_psl(output_psl, top_n)
    
    # Step 5: Save the top hits to a file
    save_top_hits(top_hits, output_hits_file)

# Example usage
input_fasta = "/lustre/scratch126/casm/team274sb/lr26/dux4_entire_gene.fasta"
reference_db = "/lustre/scratch126/casm/team274sb/lr26/T2T/T2T.fa"
output_fasta = "/lustre/scratch126/casm/team274sb/lr26/rolling_window/output_sequences.fasta"
output_psl = "/lustre/scratch126/casm/team274sb/lr26/rolling_window/blat_output.psl"
output_hits_file = "/lustre/scratch126/casm/team274sb/lr26/rolling_window/top_20_hits.txt"

process_fasta(input_fasta, reference_db, output_fasta, output_psl, output_hits_file)

