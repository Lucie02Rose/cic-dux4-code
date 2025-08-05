#!/bin/bash
### parameters fort he LSF job ###
#BSUB -n 8
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q yesterday
#BSUB -J blast-giab
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-duxblastgiab.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-duxblastgiab.err

### activate the conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define pathways to contigs ###
dir="/lustre/scratch126/casm/team274sb/lr26/Dux_search"
hg002="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG002_revio.bp.p_ctg.fasta"
hg003="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_revio.bp.p_ctg.fasta"
hg004="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG004_revio.bp.p_ctg.fasta"
hg002_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG002_combined.bp.p_ctg.fasta"
hg003_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG003_combined.bp.p_ctg.fasta"
hg004_comb="/lustre/scratch126/casm/team274sb/lr26/GIAB/HG004_combined.bp.p_ctg.fasta"

### fasta files with regions of interest ###
dux="/lustre/scratch126/casm/team274sb/lr26/Dux_search/dux4_intron2.fasta"
rpl="/lustre/scratch126/casm/team274sb/lr26/Dux_search/rpl23.fasta"
plam="/lustre/scratch126/casm/team274sb/lr26/Dux_search/dux4_plam_pa_beta_sat.fasta"
dux_rpl="/lustre/scratch126/casm/team274sb/lr26/Dux_search/dux_rpl.fasta"

### search for these in the hifi contigs ###
### sort results by bit score ###
### find where the contigs are aligned (4, 10 etc) - in IGV ###
### change to the output directory ###
cd "$dir" 

### make databases out of the contig files ###
makeblastdb -in "$hg002" -dbtype nucl -out hg002
makeblastdb -in "$hg003" -dbtype nucl -out hg003
makeblastdb -in "$hg004" -dbtype nucl -out hg004
makeblastdb -in "$hg002_comb" -dbtype nucl -out hg002_comb
makeblastdb -in "$hg003_comb" -dbtype nucl -out hg003_comb
makeblastdb -in "$hg004_comb" -dbtype nucl -out hg004_comb

### blast each fasta of interest against each database, directly sorting and saving as a txt file ###
blastn -query "$dux" -db hg002 -outfmt 6 | sort -k12,12nr > dux_i2_hg002.txt
blastn -query "$dux" -db hg003 -outfmt 6 | sort -k12,12nr > dux_i2_hg003.txt
blastn -query "$dux" -db hg004 -outfmt 6 | sort -k12,12nr > dux_i2_hg004.txt
blastn -query "$dux" -db hg002_comb -outfmt 6 | sort -k12,12nr > dux_i2_hg002_comb.txt
blastn -query "$dux" -db hg003_comb -outfmt 6 | sort -k12,12nr > dux_i2_hg003_comb.txt
blastn -query "$dux" -db hg004_comb -outfmt 6 | sort -k12,12nr > dux_i2_hg004_comb.txt

blastn -query "$rpl" -db hg002 -outfmt 6 | sort -k12,12nr > rpl_hg002.txt
blastn -query "$rpl" -db hg003 -outfmt 6 | sort -k12,12nr > rpl_hg003.txt
blastn -query "$rpl" -db hg004 -outfmt 6 | sort -k12,12nr > rpl_hg004.txt 
blastn -query "$rpl" -db hg002_comb -outfmt 6 | sort -k12,12nr > rpl_hg002_comb.txt
blastn -query "$rpl" -db hg003_comb -outfmt 6 | sort -k12,12nr > rpl_hg003_comb.txt
blastn -query "$rpl" -db hg004_comb -outfmt 6 | sort -k12,12nr > rpl_hg004_comb.txt 

blastn -query "$dux_rpl" -db hg002 -outfmt 6 | sort -k12,12nr > dux_rpl_hg002.txt
blastn -query "$dux_rpl" -db hg003 -outfmt 6 | sort -k12,12nr > dux_rpl_hg003.txt
blastn -query "$dux_rpl" -db hg004 -outfmt 6 | sort -k12,12nr > dux_rpl_hg004.txt
blastn -query "$dux_rpl" -db hg002_comb -outfmt 6 | sort -k12,12nr > dux_rpl_hg002_comb.txt
blastn -query "$dux_rpl" -db hg003_comb -outfmt 6 | sort -k12,12nr > dux_rpl_hg003_comb.txt
blastn -query "$dux_rpl" -db hg004_comb -outfmt 6 | sort -k12,12nr > dux_rpl_hg004_comb.txt

blastn -query "$plam" -db hg002 -outfmt 6 | sort -k12,12nr > plam_hg002.txt
blastn -query "$plam" -db hg003 -outfmt 6 | sort -k12,12nr > plam_hg003.txt
blastn -query "$plam" -db hg004 -outfmt 6 | sort -k12,12nr > plam_hg004.txt
blastn -query "$plam" -db hg002_comb -outfmt 6 | sort -k12,12nr > plam_hg002_comb.txt
blastn -query "$plam" -db hg003_comb -outfmt 6 | sort -k12,12nr > plam_hg003_comb.txt
blastn -query "$plam" -db hg004_comb -outfmt 6 | sort -k12,12nr > plam_hg004_comb.txt
