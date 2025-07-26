#!/bin/bash
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/wgspic.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/wgspic.e
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
#reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
output_dir="/lustre/scratch126/casm/team274sb/lr26/wgs-t2t"
input_patient_bam="/lustre/scratch126/casm/team274sb/lr26/WGS/PD54858d/PD54858d.sample.dupmarked.bam"
input_mom_bam="/lustre/scratch126/casm/team274sb/lr26/WGS/PD54859b/PD54859b.sample.dupmarked.bam"
#angus_patient="/nfs/cancer_ref01/nst_links/live/3030/PD59703a/PD59703a.sample.dupmarked.bam"
# Move to the output directory
#echo "Making output directory and moving there"
#mkdir -p "$output_dir"
cd "$output_dir"

# Extract the reads from the bam file and name them
export _JAVA_OPTIONS="-Xmx250G"

echo "Extracting paired-end and singleton reads from patient WGS"
#samtools fastq -F 0 -1 patient_reads_R1.fastq -2 patient_reads_R2.fastq -s patient_singletons.fastq "$input_patient_bam"

picard SamToFastq I="$input_patient_bam" F=R1_patient_pic.fastq F2=R2_patient_pic.fastq FU=singletons_patient_pic.fastq

picard SamToFastq I="$input_mom_bam" F=R1_mom_pic.fastq F2=R2_mom_pic.fastq FU=singletons_mom_pic_.fastq


#echo "Extracting paired-end and singleton reads from mom WGS"
#samtools fastq -F 0 -1 mom_reads_R1.fastq -2 mom_reads_R2.fastq -s mom_singletons.fastq "$input_mom_bam"

#picard SamToFastq I="$input_mom_bam" F=R1_mom_pic.fastq F2=R2_mom_pic.fastq FU=singletons_mom_pic.fastq

# Generate bwa-mem2 genome index if not there (28N GB memory, N = ref size)
#echo "Generating  bwa-mem2 genome index"
#bwa-mem2 index "$reference"

# Align paired-end reads
echo "Aligning patient paired-end reads"
#bwa-mem2 mem -t 32 "$reference" R1_patient_pic.fastq R2_patient_pic.fastq > patient_paired_pic.sam
#bwa-mem2 mem -t 32 "$reference" R1_PD59703a.fastq R2_PD59703a.fastq > PD59703a_paired.sam
#echo "Aligning mother paired-end reads"
#bwa-mem2 mem -t 32 "$reference" R1_mom_pic.fastq R2_mom_pic.fastq > mom_paired_pic.sam

# Align singleton reads
#echo "Aligning patient singleton reads"
#bwa-mem2 mem -t 32 "$reference" patient_singletons.fastq > patient_singletons.sam
#echo "Aligning mother singleton reads"
#bwa-mem2 mem -t 32 "$reference" mom_singletons.fastq > mom_singletons.sam

# Merge paired-end and singleton alignments
echo "Converting from sam to bam - patient"
#echo -e "@HD\tVN:1.6\tSO:coordinate" | cat - patient_singletons_pic.sam > patient_singletons_fixed_pic.sam
#samtools view -Sb PD59703a_paired.sam > PD59703a_paired.bam
#samtools view -Sb patient_paired_pic.sam > patient_paired_pic.bam
#samtools sort -Sb patient_paired_pic.bam > patient_paired_pic_sorted.bam
#echo "Merging patient bam files"
#samtools merge -@ 32 -f patient_merged.bam patient_paired.bam patient_singletons.bam

#echo "Converting from sam to bam - mom"
#echo -e "@HD\tVN:1.6\tSO:coordinate" | cat - mom_singletons.sam > mom_singletons_fixed.sam

#samtools view -Sb mom_paired_pic.sam > mom_paired_pic.bam
#samtools view -h -Sb mom_singletons_fixed.sam > mom_singletons.bam
#echo "Merging mom bam files"
#samtools merge -@ 32 -f mom_merged.bam mom_paired.bam mom_singletons.bam

# Sort and index the final BAM files
echo "Sorting and indexing patient bams"
#samtools sort -o patient_final_pic.bam patient_paired_pic.bam
#samtools index patient_final_pic.bam

#samtools sort -o PD59703a_final.bam PD59703a_paired.bam
#samtools index PD59703a_final.bam

#echo "Sorting and indexing patient bams"
#samtools sort -o mom_final_pic.bam mom_paired_pic.bam
#samtools index mom_final_pic.bam

# Cleanup intermediate files
#echo "Cleaning up intermediate files"
#rm patient_paired.sam patient_singletons.sam patient_paired.bam patient_singletons.bam
#rm mom_paired.sam mom_singletons.sam mom_paired.bam mom_singletons.bam


