#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q yesterday
#BSUB -J pbmm2-hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/giab-fasta-concat.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/giab-fasta-concat.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define directories ###
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
giab_dir="/lustre/scratch126/casm/team274sb/lr26/GIAB"

### change to the GIAB directory where raw bam files are ##
cd "$giab_dir"

### function to convert bam files to fasta using new samtools ###
### defining file names and running ###
convert_bam_to_fastq() {
    bamfile="$1"
    base=$(basename "$bamfile" .bam)
    fastq_out="$giab_dir/${base}.fastq"
    echo "Converting $bamfile to $fastq_out"
    samtools fastq "$bamfile" > "$fastq_out"
}
### exporting and running the function ###
export -f convert_bam_to_fastq
export giab_dir

### find the corresponding bam files in the directory and use the function in parallel ###
find "$giab_dir" -name "*HG00?-rep1.bam" | parallel -j 8 convert_bam_to_fastq
find "$giab_dir" -name "*HG00?_sprq.bam" | parallel -j 8 convert_bam_to_fastq

### combine the revio and sprq together if both present ###
### error handling for each individual ###
if [[ -f HG002-rep1.fastq && -f HG002_sprq.fastq ]]; then
    cat HG002-rep1.fastq HG002_sprq.fastq > combined_HG002.fastq
else
    echo "missing both fastqs, skipping."
fi

if [[ -f HG003-rep1.fastq && -f HG003_sprq.fastq ]]; then
    cat HG003-rep1.fastq HG003_sprq.fastq > combined_HG003.fastq
else
    echo "missing both fastqs, skipping."
fi

if [[ -f HG004-rep1.fastq && -f HG004_sprq.fastq ]]; then
    cat HG004-rep1.fastq HG004_sprq.fastq > combined_HG004.fastq
else
    echo "missing both fastqs, skipping."
fi
