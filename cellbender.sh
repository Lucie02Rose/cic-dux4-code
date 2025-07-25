#!/bin/bash
#BSUB -n 16
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q gpu-normal
#BSUB -J cellbender
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err
#BSUB -gpu "num=2:mode=shared" 
#BSUB "sleep 60; nvidia-smi"

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate cellbender

cellbender remove-background \
    --input "/lustre/scratch126/casm/team274sb/lr26/scRNA/raw_data/adata_merged_scrublet_scvi_umi_cons_doub.h5ad" \
    --output "/lustre/scratch126/casm/team274sb/lr26/scRNA/raw_data/adata_merged_scrublet_scvi_umi_cons_doub_cellbend_new.h5" \
    --cuda \
    --fpr 0.0 \
    --expected-cells 80000

echo "Cellbender completed." 
