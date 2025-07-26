#!/bin/bash
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

# Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.mmi"
giab_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio"
fastq_dir="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/fastq-revio-sprq"

#cd "$fastq_dir"

#seqtk seq -A m84039_241001_220042_s2.hifi_reads.bc2018.fastq.gz > HG002_sprq.fasta
#seqtk seq -A m84039_241002_020632_s4.hifi_reads.bc2021.fastq.gz > HG003_sprq.fasta
#seqtk seq -A m84039_241002_000337_s3.hifi_reads.bc2020.fastq.gz > HG004_sprq.fasta

cd "$giab_dir"

cat m84011_220902_175841_s1.HG002-rep1.fasta "$fastq_dir/HG002_sprq.fasta" > combined_HG002.fasta

#cat m84010_220919_232145_s1.HG004-rep1.fasta "$fastq_dir/HG003_sprq.fasta" > combined_HG003.fasta

#cat m84010_220919_235306_s2.HG003-rep1.fasta "$fastq_dir/HG004_sprq.fasta" > combined_HG004.fasta
