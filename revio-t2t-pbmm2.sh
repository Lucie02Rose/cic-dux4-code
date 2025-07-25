#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/t2talignment%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/t2talignment%J.e
#BSUB -n 18
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J pbmm2-alignment-t2t-revio
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

#!/bin/bash

# Path to raw FASTA reference (instead of prebuilt .mmi)
reference_fa="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
reference_mmi="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"

input_dir="/lustre/scratch126/casm/team274sb/lr26/Revio_raw"
output_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-t2t-hifi"

mkdir -p "$output_dir"
tmp_dir="$output_dir/tmp"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

# Check if .mmi index exists; if not, create it
if [ ! -f "$reference_mmi" ]; then
    echo "Creating minimap2 index for reference FASTA..."
    pbmm2 index "$reference_fa" "$reference_mmi"
    echo "Index created at $reference_mmi"
else
    echo "Reference index already exists at $reference_mmi"
fi

threads=18
parallel_jobs=3

align_and_index() {
    bam_file="$1"
    base_name=$(basename "$bam_file" .bam)
    output_bam="$output_dir/${base_name}_pbmm2-farm22-bam.bam"

    echo "Aligning $bam_file to $reference_mmi..."
    pbmm2 align "$reference_mmi" "$bam_file" "$output_bam" --preset HIFI --sort -j "$threads"

    echo "Completed $output_bam"
}

export -f align_and_index
export reference_mmi output_dir threads

find "$input_dir" -maxdepth 1 -type f -name "*hifi_reads*.bam" | \
    xargs -n 1 -P "$parallel_jobs" -I {} bash -c 'align_and_index "$@"' _ {}

echo "All alignments completed."
