#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q yesterday
#BSUB -J pbmm2-hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/giab-fasta.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/giab-fasta.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
giab_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio"

cd "$giab_dir"

NUM_THREADS=8

# Export function to run in parallel
convert_bam_to_fasta() {
    bamfile="$1"
    base=$(basename "$bamfile" .bam)
    fasta_out="$giab_dir/${base}.fasta"

    echo "Converting $bamfile → $fasta_out"
    samtools fasta "$bamfile" > "$fasta_out"
}

export -f convert_bam_to_fasta
export giab_dir

# Find matching files and run in parallel
find "$giab_dir" -name "*HG002-rep?.bam" | parallel -j "$NUM_THREADS" convert_bam_to_fasta



