#!/bin/bash
### parameters fort he LSF job ###
#BSUB -n 8
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q yesterday
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-duxblast.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-duxblast.err

### activate the conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define pathways to contigs ###
dir="/lustre/scratch126/cellgen/behjati/lr26/Dux_search"
mom="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/mom.bp.p_ctg.fasta"
blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/blood.bp.p_ctg.fasta"
revio="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/PacBio-revio.bp.p_ctg.fasta"

### fasta files with regions of interest ###
dux="/lustre/scratch126/cellgen/behjati/lr26/Dux_search/dux4_intron2.fasta"
rpl="/lustre/scratch126/cellgen/behjati/lr26/Dux_search/rpl23.fasta"
plam="/lustre/scratch126/cellgen/behjati/lr26/Dux_search/dux4_plam_pa_beta_sat.fasta"
dux_rpl="/lustre/scratch126/cellgen/behjati/lr26/Dux_search/dux_rpl.fasta"

### search for these in the hifi contigs ###
### sort results by bit score ###
### find where the contigs are aligned (4, 10 etc) - in IGV ###
### change to the output directory ###
cd "$dir" 

### make databases out of the contig files ###
makeblastdb -in "$mom" -dbtype nucl -out hifi_mom
makeblastdb -in "$blood" -dbtype nucl -out hifi_blood
makeblastdb -in "$tumor" -dbtype nucl -out hifi_tumor

### blast each fasta of interest against each database, directly sorting and saving as a txt file ###
blastn -query "$dux" -db hifi_mom -outfmt 6 | sort -k12,12nr > dux_i2_hifimom.txt
blastn -query "$dux" -db hifi_blood -outfmt 6 | sort -k12,12nr > dux_i2_hifiblood.txt
blastn -query "$dux" -db hifi_tumor -outfmt 6 | sort -k12,12nr > dux_i2_hifitum.txt

blastn -query "$rpl" -db hifi_mom -outfmt 6 | sort -k12,12nr > rpl_hifimom.txt
blastn -query "$rpl" -db hifi_blood -outfmt 6 | sort -k12,12nr > rpl_hifiblood.txt
blastn -query "$rpl" -db hifi_tumor -outfmt 6 | sort -k12,12nr > rpl_hifitum.txt 

blastn -query "$dux_rpl" -db hifi_mom -outfmt 6 | sort -k12,12nr > dux_rpl_hifimom.txt
blastn -query "$dux_rpl" -db hifi_blood -outfmt 6 | sort -k12,12nr > dux_rpl_hifiblood.txt
blastn -query "$dux_rpl" -db hifi_tumor -outfmt 6 | sort -k12,12nr > dux_rpl_hifitum.txt

blastn -query "$plam" -db hifi_mom -outfmt 6 | sort -k12,12nr > plam_hifimom.txt
blastn -query "$plam" -db hifi_blood -outfmt 6 | sort -k12,12nr > plam_hifiblood.txt
blastn -query "$plam" -db hifi_tumor -outfmt 6 | sort -k12,12nr > plam_hifitum.txt


