#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q normal
#BSUB -J rna-blats
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-rna-blat.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-rna-blat.e
### note that this script was used and written before the farmageddon ###
### activate the conda base environment (not really necessary) ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### directories to use ###
reference="/nfs/users/nfs_l/lr26/PacBio/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"
output_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t/blats"
### executable blat ###
blat="/nfs/users/nfs_l/lr26/blat"

### make the output directory and move there ### 
mkdir -p "$output_dir"
cd "$output_dir"

### fasta files for the CIC-DUX4 cohort patients
fastas=("GSM5024895" "GSM5024897" "GSM5024898")
### blat function for each fasta file and sequence ###
blat_exec() {
	fasta=$1
	fasta_file="${input_dir}/${fasta}.fasta"
	psl_out="${output_dir}/${fasta}.psl"
	if [[ -f "$fasta_file" ]]; then
		echo "Processing $fasta"
		"$blat" "$reference" "$fasta_file" "$psl_out" -out=psl -tileSize=7 -minMatch=1 -oneOff=1 -repMatch=1000000 -maxIntron=1000000 -stepSize=3 -minScore=10 
		echo "Processed $fasta"
	else
		echo "Problem - no fasta found"
	fi
}
### do for all the fasta files ###
for fasta in "${fastas[@]}"; do
	blat_exec "$fasta"
done



