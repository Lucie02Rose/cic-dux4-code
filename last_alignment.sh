#!/bin/bash
#BSUB -n 16
#BSUB -M 70000
#BSUB -R 'span[hosts=1] select[mem>70000] rusage[mem=70000]'
#BSUB -q normal
#BSUB -J last
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-last.o
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-last.e

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories
input_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-last"
reference="/lustre/scratch126/cellgen/behjati/lr26/last/"
db_prefix="t2t_db"
input_file="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/blood_1C01_hifi_reads.fastq.gz"

# List of gzipped FASTQ files
samples=("tumor_1B01_hifi_reads.fastq.gz" "mom_1B02_hifi_reads.fastq.gz" "blood_1C01_hifi_reads.fastq.gz" "tumor_all_hifi_reads.fastq.gz")

# Train model once — make sure $train_input is set to one FASTQ.gz to convert for training
train_input="$input_dir/${samples[0]}"
train_fasta="$reference/train_input.fasta"

# Convert training input FASTQ.gz to FASTA for last-train
zcat "$train_input" | seqtk seq -a - > "$train_fasta"
last-train -P16 "$reference/$db_prefix" "$train_fasta" > "$reference/model.train"

# Loop over each sample
for sample in "${samples[@]}"; do
  echo "Processing $sample"

  input_fastq="$input_dir/$sample"
  base_name="${sample%.fastq.gz}"
  output_base="$output_dir/$base_name"
  sample_fasta="$output_dir/${base_name}.fasta"

  # Convert FASTQ.gz to FASTA for this sample
  zcat "$input_fastq" | seqtk seq -a - > "$sample_fasta"

  # Run LAST alignment on the FASTA
  lastal -P16 -C2 -p "$reference/model.train" "$reference/$db_prefix" "$sample_fasta" | last-split > "${output_base}.maf"

  # Optional: convert MAF to sorted BAM for IGV
  maf-convert sam "${output_base}.maf" | samtools view -bS - | samtools sort -o "${output_base}.sorted.bam"
  samtools index "${output_base}.sorted.bam"
done
