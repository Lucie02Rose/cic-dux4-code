#!/bin/bash
#BSUB -o /lustre/scratch125/cellgen/behjati/lr26/outputs/%J-wgspic.o
#BSUB -e /lustre/scratch125/cellgen/behjati/lr26/errors/%J-wgspic.e
#BSUB -n 32
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q hugemem
#BSUB -J wgs-t2t-alignment-patient
#BSUB -G team274

### Activate conda environment
echo "Activating bwa-mem2 environment"
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bwa-mem2

### Directories to process
reference="/lustre/scratch125/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
output_dir="/lustre/scratch125/cellgen/behjati/lr26/WGS"
fastqc="/lustre/scratch125/cellgen/behjati/lr26/WGS/fastqc"
input_tumor_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54858d.v1.sample.dupmarked.bam"
input_mom_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54859b.v1.sample.dupmarked.bam"
input_blood_bam="/lustre/scratch125/cellgen/behjati/lr26/WGS/PD54858b.v1.sample.dupmarked.bam"
#angus_patient="/nfs/cancer_ref01/nst_links/live/3030/PD59703a/PD59703a.sample.dupmarked.bam"
# Move to the output directory
#echo "Making output directory and moving there"
#mkdir -p "$output_dir"
cd "$output_dir"
mkdir -p "$fastqc"

# Extract the reads from the bam file and name them
export _JAVA_OPTIONS="-Xmx250G"

#echo "Extracting paired-end and singleton reads from patient WGS"

#picard SamToFastq I="$input_tumor_bam" F=R1_tumor_pic.fastq F2=R2_tumor_pic.fastq FU=singletons_tumor_pic.fastq

#picard SamToFastq I="$input_blood_bam" F=R1_blood_pic.fastq F2=R2_blood_pic.fastq FU=singletons_blood_pic.fastq

#picard SamToFastq I="$input_mom_bam" F=R1_mom_pic.fastq F2=R2_mom_pic.fastq FU=singletons_mom_pic.fastq

#conda deactivate
#conda activate base

#find "$output_dir" -name "*.fastq" | \
 #   parallel -j 8 fastqc {} -o "$fastqc"

# Generate bwa-mem2 genome index if not there (28N GB memory, N = ref size)
#conda deactivate 
#conda activate bwa-mem2
echo "Generating  bwa-mem2 genome index"
bwa-mem2 index "$reference"

# Align paired-end reads
echo "Aligning patient paired-end reads"
bwa-mem2 mem -t 32 "$reference" R1_blood_pic.fastq R2_blood_pic.fastq > blood_paired_pic.sam
bwa-mem2 mem -t 32 "$reference" R1_tumor_pic.fastq R2_tumor_pic.fastq > tumor_paired_pic.sam
#echo "Aligning mother paired-end reads"
bwa-mem2 mem -t 32 "$reference" R1_mom_pic.fastq R2_mom_pic.fastq > mom_paired_pic.sam

# Align singleton reads
echo "Aligning patient singleton reads"
bwa-mem2 mem -t 32 "$reference" tumor_singletons.fastq > tumor_singletons.sam
bwa-mem2 mem -t 32 "$reference" blood_singletons.fastq > blood_singletons.sam
echo "Aligning mother singleton reads"
bwa-mem2 mem -t 32 "$reference" mom_singletons.fastq > mom_singletons.sam

# Merge paired-end and singleton alignments
echo "Converting from sam to bam - patient"
echo -e "@HD\tVN:1.6\tSO:coordinate" | cat - blood_singletons_pic.sam > blood_singletons_fixed_pic.sam
samtools view -Sb blood_paired_pic.sam > blood_paired_pic.bam
samtools sort -Sb blood_paired_pic.bam > blood_paired_pic_sorted.bam
echo "Merging patient bam files"
samtools merge -@ 32 -f blood_merged.bam blood_paired.bam blood_singletons.bam

echo "Converting from sam to bam - patient"
echo -e "@HD\tVN:1.6\tSO:coordinate" | cat - tumor_singletons_pic.sam > tumor_singletons_fixed_pic.sam
samtools view -Sb tumor_paired_pic.sam > tumor_paired_pic.bam
samtools sort -Sb tumor_paired_pic.bam > tumor_paired_pic_sorted.bam
echo "Merging patient bam files"
samtools merge -@ 32 -f tumor_merged.bam tumor_paired.bam tumor_singletons.bam

echo "Converting from sam to bam - mom"
echo -e "@HD\tVN:1.6\tSO:coordinate" | cat - mom_singletons.sam > mom_singletons_fixed.sam
samtools view -Sb mom_paired_pic.sam > mom_paired_pic.bam
samtools view -h -Sb mom_singletons_fixed.sam > mom_singletons.bam
echo "Merging mom bam files"
samtools merge -@ 32 -f mom_merged.bam mom_paired.bam mom_singletons.bam

# Sort and index the final BAM files
echo "Sorting and indexing patient bams"
samtools sort -o blood_final_pic.bam blood_paired_pic.bam
samtools index blood_final_pic.bam
samtools sort -o tumor_final_pic.bam tumor_paired_pic.bam
samtools index tumor_final_pic.bam

echo "Sorting and indexing mom bams"
samtools sort -o mom_final_pic.bam mom_paired_pic.bam
samtools index mom_final_pic.bam

# Cleanup intermediate files
#echo "Cleaning up intermediate files"
#rm patient_paired.sam patient_singletons.sam patient_paired.bam patient_singletons.bam
#rm mom_paired.sam mom_singletons.sam mom_paired.bam mom_singletons.bam


