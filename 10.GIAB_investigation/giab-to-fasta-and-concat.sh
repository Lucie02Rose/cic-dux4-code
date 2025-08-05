#!/bin/bash
### parameters for the lsf job ###
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q yesterday
#BSUB -J pbmm2-hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/giab-fasta-concat.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/giab-fasta-concat.err

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
convert_bam_to_fasta() {
    bamfile="$1"
    base=$(basename "$bamfile" .bam)
    fasta_out="$giab_dir/${base}.fasta"
    echo "Converting $bamfile to $fasta_out"
    samtools fasta "$bamfile" > "$fasta_out"
}
### exporting and running the function ###
export -f convert_bam_to_fasta
export giab_dir

### find the corresponding bam files in the directory and use the function in parallel ###
find "$giab_dir" -name "*HG00?-rep?.bam" | parallel -j 8 convert_bam_to_fasta
find "$giab_dir" -name "*HG00?_sprq.bam" | parallel -j 8 convert_bam_to_fasta

### combine the revio and sprq together if both present ###
### error handling for each individual ###
if [[ -f HG002-rep1.fasta && -f HG002_sprq.fasta ]]; then
    cat HG002-rep1.fasta HG002_sprq.fasta > combined_HG002.fasta
else
    echo "WARNING: Missing FASTA(s) for HG002. Skipping concatenation."
fi

if [[ -f HG003-rep1.fasta && -f HG003_sprq.fasta ]]; then
    cat HG003-rep1.fasta HG003_sprq.fasta > combined_HG003.fasta
else
    echo "WARNING: Missing FASTA(s) for HG003. Skipping concatenation."
fi

if [[ -f HG004-rep1.fasta && -f HG004_sprq.fasta ]]; then
    cat HG004-rep1.fasta HG004_sprq.fasta > combined_HG004.fasta
else
    echo "WARNING: Missing FASTA(s) for HG004. Skipping concatenation."
fi
