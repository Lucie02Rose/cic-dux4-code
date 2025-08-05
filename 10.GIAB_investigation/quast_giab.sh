#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q normal
#BSUB -J quast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-quast-giab.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-quast.err

### activate the quast conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate quast_busco

### define directories fo all files to evaluate ###
hg002="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG002_revio.bp.p_ctg.fasta"
hg003="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_revio.bp.p_ctg.fasta"
hg004="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG004_revio.bp.p_ctg.fasta"
hg002_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG002_combined.bp.p_ctg.fasta"
hg003_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_combined.bp.p_ctg.fasta"
hg004_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG004_combined.bp.p_ctg.fasta"

### define and make the output directory ###
output_dir="/lustre/scratch126/casm/team274sb/lr26/GIAB/Quast"

mkdir -p "$output_dir"

### run quast on all assemblies ###
quast -o "$output_dir" "$hg002" "$hg003" "$hg004" "$hg002_comb" "$hg003_comb" "$hg004_comb"
