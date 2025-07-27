#!/bin/bash
### prameters for the LSF job ###
#BSUB -n 32
#BSUB -M 250000
#BSUB -R 'span[hosts=1] select[mem>250000] rusage[mem=250000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-aligtumor.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-aligtumor.err

### activate the base conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### list of files - reference, fastq, combined and input output directories ###
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.mmi"
input_fastq1="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A01_hifi_reads.fastq.gz"
input_fastq2="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1A02_hifi_reads.fastq.gz"
input_fastq3="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_2B01_hifi_reads.fastq.gz"
input_fastq4="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_all_4_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/tmp"

### combine all the fastq files ###
cat "$input_fastq1" "$input_fastq2" "$input_fastq3" "$input_fastq4" > "$combined_fastq"

### create the temporary and output directories if not already present ###
mkdir -p "$tmp_dir" "$output_dir"

### define the output bam file name ###
base_name=$(basename "$combined_fastq" .fastq.gz)
output_bam="$output_dir/${base_name}_pbmm2.bam"

#### temporary directory for sorting (export) ###
export TMPDIR="$tmp_dir"

### run the pbmm2 alignmment on the fastq file ###
echo "Aligning $input_fastq to $reference..."
pbmm2 align "$reference" "$combined_fastq" "$output_bam" --preset HIFI --sort -j 32 --unmapped

### end process message ###
echo "pbmm2 alignment completed: $output_bam"
