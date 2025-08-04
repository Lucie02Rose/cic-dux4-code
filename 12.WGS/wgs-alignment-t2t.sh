#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 32
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q hugemem
#BSUB -J wgs-t2t-alignment-patient
#BSUB -G team274
#BSUB -o /lustre/scratch125/cellgen/behjati/lr26/outputs/%J-wgspic.o
#BSUB -e /lustre/scratch125/cellgen/behjati/lr26/errors/%J-wgspic.e

### activate the bwa-mem2 environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bwa-mem2

### which directories to process ###
reference="/lustre/scratch125/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
output_dir="/lustre/scratch125/cellgen/behjati/lr26/WGS"
fastqc="/lustre/scratch125/cellgen/behjati/lr26/WGS/fastqc"
input_tumor_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54858d.v1.sample.dupmarked.bam"
input_mom_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54859b.v1.sample.dupmarked.bam"
input_blood_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54858b.v1.sample.dupmarked.bam"

### make and move to directories ###
echo "Making output directory and moving there"
mkdir -p "$output_dir"
cd "$output_dir"
mkdir -p "$fastqc"

### export the java memory needed ###
export _JAVA_OPTIONS="-Xmx250G"
### extract using picard ###
echo "Extracting paired-end and singleton reads from patient WGS"
picard SamToFastq I="$input_tumor_bam" F=R1_tumor_pic.fastq F2=R2_tumor_pic.fastq FU=singletons_tumor_pic.fastq
picard SamToFastq I="$input_blood_bam" F=R1_blood_pic.fastq F2=R2_blood_pic.fastq FU=singletons_blood_pic.fastq
picard SamToFastq I="$input_mom_bam" F=R1_mom_pic.fastq F2=R2_mom_pic.fastq FU=singletons_mom_pic.fastq
### deactivate and activate base again ###
conda deactivate
conda activate base
### run fastqc ###
find "$output_dir" -name "*.fastq" | \
    parallel -j 8 fastqc {} -o "$fastqc"
    
### activate bwa-mem2 and generate the reference index ###
conda deactivate 
conda activate bwa-mem2

echo "Generating  bwa-mem2 genome index"
bwa-mem2 index "$reference"

### align paired-end reads ###
echo "Aligning patient paired-end reads"
bwa-mem2 mem -t 32 "$reference" R1_blood_pic.fastq R2_blood_pic.fastq > blood_paired_pic.sam
bwa-mem2 mem -t 32 "$reference" R1_tumor_pic.fastq R2_tumor_pic.fastq > tumor_paired_pic.sam
bwa-mem2 mem -t 32 "$reference" R1_mom_pic.fastq R2_mom_pic.fastq > mom_paired_pic.sam

### process blood, tumor mom reads - sort and index by samtools ###
echo "Processing blood"
samtools sort -o blood_final_pic.bam blood_paired_pic.sam
samtools index blood_final_pic.bam
echo "Processing tumor"
samtools sort -o tumor_final_pic.bam tumor_paired_pic.sam
samtools index tumor_final_pic.bam
echo "Processing mom"
samtools sort -o mom_final_pic.bam mom_paired_pic.sam
samtools index mom_final_pic.bam

### clean up all the files that are not needed ###
echo "Cleaning up intermediate files"
rm blood_paired_pic.sam tumor_paired_pic.sam mom_paired_pic.sam
### these are usually empty ###
rm singletons_tumor_pic.fastq singletons_blood_pic.fastq singletons_mom_pic.fastq
