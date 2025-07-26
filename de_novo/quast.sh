#!/bin/bash
#BSUB -n 16
#BSUB -M 50000
#BSUB -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'
#BSUB -q normal
#BSUB -J quast
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/quast%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/quast%J.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate quast_busco

### Directories
mom_hap1="/lustre/scratch126/casm/team274sb/lr26/hifiasm_mom_denovo/mom_denovo_hifiasm.bp.hap1.p_ctg.fasta"
mom_hap2="/lustre/scratch126/casm/team274sb/lr26/hifiasm_mom_denovo/mom_denovo_hifiasm.bp.hap2.p_ctg.fasta"
mom="/lustre/scratch126/casm/team274sb/lr26/hifiasm_mom_denovo/mom_denovo_hifiasm.bp.p_ctg.fasta"

blood_hap1="/lustre/scratch126/casm/team274sb/lr26/hifiasm_blood_denovo/blood_denovo_hifiasm.bp.hap1.p_ctg.fasta"
blood_hap2="/lustre/scratch126/casm/team274sb/lr26/hifiasm_blood_denovo/blood_denovo_hifiasm.bp.hap2.p_ctg.fasta"
blood="/lustre/scratch126/casm/team274sb/lr26/hifiasm_blood_denovo/blood_denovo_hifiasm.bp.p_ctg.fasta"

tumor_hap1="/lustre/scratch126/casm/team274sb/lr26/hifiasm_tumor_denovo/tumor_denovo_hifiasm.bp.hap1.p_ctg.fasta"
tumor_hap2="/lustre/scratch126/casm/team274sb/lr26/hifiasm_tumor_denovo/tumor_denovo_hifiasm.bp.hap2.p_ctg.fasta"
tumor="/lustre/scratch126/casm/team274sb/lr26/hifiasm_tumor_denovo/tumor_denovo_hifiasm.bp.p_ctg.fasta"

hg002_hap1="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG002_denovo_hifiasm.bp.hap1.p_ctg.fasta"
hg002_hap2="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG002_denovo_hifiasm.bp.hap2.p_ctg.fasta"
hg002="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG002_denovo_hifiasm.bp.p_ctg.fasta"

hg003_hap1="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_denovo_hifiasm.bp.hap1.p_ctg.fasta"
hg003_hap2="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_denovo_hifiasm.bp.hap2.p_ctg.fasta"
hg003="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG003_denovo_hifiasm.bp.p_ctg.fasta"

hg004_hap1="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_denovo_hifiasm.bp.hap1.p_ctg.fasta"
hg004_hap2="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_denovo_hifiasm.bp.hap2.p_ctg.fasta"
hg004="/lustre/scratch126/casm/team274sb/lr26/pbmm2-giab-revio/HG004_denovo_hifiasm.bp.p_ctg.fasta"

# Define the output directory
output_dir="/lustre/scratch126/casm/team274sb/lr26/hifiasm_quast_output_nonsprq"

#mkdir -p "$output_dir"

# Run QUAST on all assemblies
quast -o "$output_dir" "$hg002_hap1" "$hg002_hap2" "$hg002" "$hg003_hap1" "$hg003_hap2" "$hg003" "$hg004_hap1" "$hg004_hap2" "$hg004"






