#!/bin/bash
#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J ngmlr-alignment-t2t-lucie-1A01-new
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_fastq="/lustre/scratch126/casm/team274sb/lr26/fastq/m64094e_230126_154129.hifi_reads.fastq.fastq.gz"
output_bam="/lustre/scratch126/casm/team274sb/lr26/ngmlr-alignment/1_A01/m64094e_230126_154129_ngmlr-new-mult.bam"

### Run NGMLR alignment

echo "Aligning $input_bam to $reference..."
ngmlr -r "$reference" -q "$input_fastq" -o "$output_bam" -t 16 -x pacbio

echo "NGMLR alignment completed: $output_bam"
