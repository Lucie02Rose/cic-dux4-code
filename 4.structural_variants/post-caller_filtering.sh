#!/bin/bash
### parameters fort eh LSF job ###
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-svfilt.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-svfilt.e
#BSUB -n 12
#BSUB -M 32000
#BSUB -R 'span [hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J svfilt
#BSUB -G team274

### activate the conda environment with bcftools ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### directories to process ###
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish"
tumor_all="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/tumor_all_4_hifi_reads_pbmm2_joint_call/genotyped.sv.vcf.gz"
blood="/lustre/scratch126/cellgen/behjati/lr26/PacBio-sawfish/blood_1C01_hifi_reads_pbmm2_joint_call/genotyped.sv.vcf.gz"
### note that this is not the nanomonsv output directory but I have created the PacBio-nanomonsv directory and moved my somatic tumor file here prior for consistency among variant callers ###
nanomon_un="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.vcf" 
nanomon="/lustre/scratch126/cellgen/behjati/lr26/PacBio-nanomonsv/tumor-all.nanomonsv.result.sorted.vcf"

### the reason for this script is to make variant caller ouputs fit for annotation ###
### sawfish needs filtering to retain only the somatic variants in the tumor ###
### change to the output directory ###
cd "$output_dir"
### normalise the sawfish output ###
bcftools norm -m -any -Oz -o tumor_all_genotyped.norm.sv.vcf.gz "$tumor_all"
bcftools norm -m -any -Oz -o blood_genotyped.norm.sv.vcf.gz "$blood"
### sort the sawfish output ###
bcftools sort tumor_all_genotyped.norm.sv.vcf.gz -Oz -o tumor_all_genotyped.norm.sort.sv.vcf.gz
bcftools sort blood_genotyped.norm.sv.vcf.gz -Oz -o blood_genotyped.norm.sort.sv.vcf.gz
### index the sawfish output ###
bcftools index tumor_all_genotyped.norm.sort.sv.vcf.gz
bcftools index blood_genotyped.norm.sort.sv.vcf.gz
### filter what is unique to the tumor (keep somatic variants ###
bcftools isec -p isec_tumor_only -C tumor_all_genotyped.norm.sort.sv.vcf.gz blood_genotyped.norm.sort.sv.vcf.gz

### nanomonsv output is not sorted, meaning that I have to sort it before doing anything with the file ###
bcftools sort "$nanomon_un" -o "$nanomon"

### normalisation, sorting and indexing should fix any downstream issues ###
### however if there are errors e.g. malformed starts and end positions (actually might be an inversion) ###
### try querying the ones where the start iw bigger than the end, e.g. awk '$3 < $2 {print $1"\t"$2}' ###
### label those as bad variants in a txt file and then using view exclude all the ones from the file cleaning the file ###

