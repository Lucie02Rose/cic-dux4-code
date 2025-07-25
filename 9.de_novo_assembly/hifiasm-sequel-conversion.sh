#!/bin/bash
#BSUB -n 16
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q yesterday
#BSUB -J pbmm2-hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-giabhificon.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-giabhificon.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

# Directories
reference="/nfs/users/nfs_l/lr26/nextflow_pipeline/reference/T2T/chm13v2.0.mmi"
giab_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-denovo"

cd "$giab_dir"

for file in *.p_ctg.gfa; do
    awk '/^S/{print ">"$2;print $3}' "$file" > "${file%.gfa}.fasta"
done

