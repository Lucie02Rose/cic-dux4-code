#!/bin/bash
### parameters fort he LSF job ###
#BSUB -n 8
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q yesterday
#BSUB -J blast
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-duxblastgiab.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-duxblastgiab.err

### activate the conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Define paths
dir="/lustre/scratch126/cellgen/behjati/lr26/Dux_search"
mom="/lustre/scratch126/cellgen/behjati/lr26/PacBio-mom/mom.bp.p_ctg.fasta"
blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-blood/blood_denovo_hifiasm.bp.p_ctg.fasta"
revio="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/tumor_denovo_hifiasm.bp.p_ctg.fasta"

hg002="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG002_new.bp.p_ctg.fasta"
hg003="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_new.bp.p_ctg.fasta"
hg004="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_new.bp.p_ctg.fasta"

dux="/lustre/scratch126/casm/team274sb/lr26/dux_search/dux4_intron2.fasta"
rpl="/lustre/scratch126/casm/team274sb/lr26/dux_search/rpl23.fasta"
plam="/lustre/scratch126/casm/team274sb/lr26/dux_search/dux4_plam_pa_beta_sat.fasta"
dux_rpl="/lustre/scratch126/casm/team274sb/lr26/dux_search/dux_rpl.fasta"

dux50="/lustre/scratch126/casm/team274sb/lr26/dux_search/dux4_intron2_window50.fasta"
rpl50="/lustre/scratch126/casm/team274sb/lr26/dux_search/rpl23_window50.fasta"

dux100="/lustre/scratch126/casm/team274sb/lr26/dux_search/dux4_window100.fasta"
rpl100="/lustre/scratch126/casm/team274sb/lr26/dux_search/rpl23_window100.fasta"

# search for these two in the hifi contigs
# sort by bit score

cd "$dir" 

#makeblastdb -in "$hg002" -dbtype nucl -out hifi_hg002_new
#makeblastdb -in "$hg003" -dbtype nucl -out hifi_hg003_new
#makeblastdb -in "$hg004" -dbtype nucl -out hifi_hg004_new

#makeblastdb -in "$mom" -dbtype nucl -out hifi_mom
#makeblastdb -in "$blood" -dbtype nucl -out hifi_blood
#makeblastdb -in "$tumor" -dbtype nucl -out hifi_tumor

#blastn -query "$dux" -db hifi_hg002_new -outfmt 6 | sort -k12,12nr > dux_i2_hifihg002_new.txt

#blastn -query "$dux" -db hifi_hg003_new -outfmt 6 | sort -k12,12nr > dux_i2_hifihg003_new.txt
#blastn -query "$dux" -db hifi_hg004_new -outfmt 6 | sort -k12,12nr > dux_i2_hifihg004_new.txt

#blastn -query "$dux" -db hifi_mom -outfmt 6 | sort -k12,12nr > dux_i2_hifimom.txt
#blastn -query "$dux" -db hifi_blood -outfmt 6 | sort -k12,12nr > dux_i2_hifiblood.txt
#blastn -query "$dux" -db hifi_tumor -outfmt 6 | sort -k12,12nr > dux_i2_hifitum.txt

#blastn -query "$rpl" -db hifi_hg002_new -outfmt 6 | sort -k12,12nr > rpl_hifihg002_new.txt
#blastn -query "$rpl" -db hifi_hg003_new -outfmt 6 | sort -k12,12nr > rpl_hifihg003_new.txt
#blastn -query "$rpl" -db hifi_hg004_new -outfmt 6 | sort -k12,12nr > rpl_hifihg004_new.txt

#blastn -query "$rpl" -db hifi_mom -outfmt 6 | sort -k12,12nr > rpl_hifimom.txt
#blastn -query "$rpl" -db hifi_blood -outfmt 6 | sort -k12,12nr > rpl_hifiblood.txt
#blastn -query "$rpl" -db hifi_tumor -outfmt 6 | sort -k12,12nr > rpl_hifitum.txt 

#blastn -query "$dux_rpl" -db hifi_hg002_new -outfmt 6 | sort -k12,12nr > dux_rpl_hifihg002_new.txt

#blastn -query "$dux_rpl" -db hifi_hg003_new -outfmt 6 | sort -k12,12nr > dux_rpl_hifihg003_new.txt
#blastn -query "$dux_rpl" -db hifi_hg004_new -outfmt 6 | sort -k12,12nr > dux_rpl_hifihg004_new.txt

#blastn -query "$dux_rpl" -db hifi_mom -outfmt 6 | sort -k12,12nr > dux_rpl_hifimom.txt
#blastn -query "$dux_rpl" -db hifi_blood -outfmt 6 | sort -k12,12nr > dux_rpl_hifiblood.txt
#blastn -query "$dux_rpl" -db hifi_tumor -outfmt 6 | sort -k12,12nr > dux_rpl_hifitum.txt

blastn -query "$plam" -db hifi_hg002_new -outfmt 6 | sort -k12,12nr > plam_hifihg002_new.txt

blastn -query "$plam" -db hifi_hg003_new -outfmt 6 | sort -k12,12nr > plam_hifihg003_new.txt
blastn -query "$plam" -db hifi_hg004_new -outfmt 6 | sort -k12,12nr > plam_hifihg004_new.txt

#blastn -query "$plam" -db hifi_mom -outfmt 6 | sort -k12,12nr > plam_hifimom.txt
#blastn -query "$plam" -db hifi_blood -outfmt 6 | sort -k12,12nr > plam_hifiblood.txt
#blastn -query "$plam" -db hifi_tumor -outfmt 6 | sort -k12,12nr > plam_hifitum.txt

# then search for 50bp windows in the RNA data and 100bp windows in the WGS data
#seqtk seq -A our_patient_final_merged_R1.fastq > RNA_patient_R1.fasta
#seqtk seq -A our_patient_final_merged_R2.fastq > RNA_patient_R2.fasta

#cat RNA_patient_R1.fasta RNA_patient_R2.fasta > RNA_patient_combined.fasta

#makeblastdb -in RNA_patient_combined.fasta -dbtype nucl -out rna_patient_db

#blastn -query "$dux50" -db rna_patient_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=1 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > dux_rnatum_top1.txt

#blastn -query "$rpl50" -db rna_patient_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=1 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > rpl_rnatum_top1.txt

#seqtk seq -A /lustre/scratch126/casm/team274sb/lr26/wgs-t2t/R1_patient_pic.fastq > DNA_patient_R1.fasta
#seqtk seq -A /lustre/scratch126/casm/team274sb/lr26/wgs-t2t/R2_patient_pic.fastq > DNA_patient_R2.fasta
#seqtk seq -A /lustre/scratch126/casm/team274sb/lr26/wgs-t2t/R1_mom_pic.fastq > DNA_mom_R1.fasta
#seqtk seq -A /lustre/scratch126/casm/team274sb/lr26/wgs-t2t/R2_mom_pic.fastq > DNA_mom_R2.fasta

#cat DNA_patient_R1.fasta DNA_patient_R2.fasta > DNA_patient_combined.fasta
#cat DNA_mom_R1.fasta DNA_mom_R2.fasta > DNA_mom_combined.fasta

#makeblastdb -in DNA_patient_combined.fasta -dbtype nucl -out dna_patient_db
#makeblastdb -in DNA_mom_combined.fasta -dbtype nucl -out dna_mom_db

#blastn -query "$dux100" -db dna_patient_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=2 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > dux_dnablood_top2.txt

#blastn -query "$rpl100" -db dna_patient_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=2 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > rpl_dnablood_top2.txt

#blastn -query "$dux100" -db dna_mom_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=2 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > dux_dnamom_top2.txt

#blastn -query "$rpl100" -db dna_mom_db -outfmt 6 -max_hsps 1 -num_threads 4 |
#sort -k1,1 -k12,12nr | awk -v max=2 '{
#  count[$1]++;
#  if (count[$1] <= max) print;
#}' > rpl_dnamom_top2.txt


