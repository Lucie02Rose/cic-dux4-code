#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q long
#BSUB -J improv_allele-counter
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
positions_bed="/lustre/scratch126/casm/team274sb/lr26/allele-integrator-pbfix/GenotypingResults/patient_output.bed"
s1A01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A01/m64094e_230126_154129.hifi_reads_pbmm2-farm22-bam.bam"
s1A02="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A02/m64178e_230206_134948.hifi_reads_pbmm2-farm22-bam.bam"
s2B01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/2_B01/m64178e_230207_165902.hifi_reads_pbmm2-farm22-bam.bam"
s1B01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B01/m84047_230404_172053_s2.hifi_reads.default_pbmm2-farm22-bam.bam"
s1B02="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B02/m84047_240202_152510_s2.hifi_reads.bc2025_pbmm2-farm22-bam.bam"
s1C01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_C01/m84047_240202_155616_s3.hifi_reads.bc2026_pbmm2-farm22-bam.bam"
mom="/nfs/cancer_ref01/nst_links/live/3306/PD54859b/PD54859b.sample.dupmarked.bam"
patient="/nfs/cancer_ref01/nst_links/live/3306/PD54858d/PD54858d.sample.dupmarked.bam"

set -exo pipefail

samtools mpileup -d 200 -f "$reference" -l "$positions_bed" "$s1C01" > s1C01_raw.tsv

samtools mpileup -f "$reference" -l "$positions_bed" "$patient" > patient_wgs_raw.tsv

samtools mpileup -f "$reference" -l "$positions_bed" "$mom" > mom_wgs_raw.tsv














