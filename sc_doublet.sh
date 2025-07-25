#!/bin/bash
#BSUB -n 8
#BSUB -M 64000
#BSUB -R 'span[hosts=1] select[mem>64000] rusage[mem=64000]'
#BSUB -q normal
#BSUB -J sc-doublets
#BSUB -G team274
#BSUB -o /nfs/users/nfs_l/lr26/errors/%J-sc-doublets.out
#BSUB -e /nfs/users/nfs_l/lr26/errors/%J-sc-doublets.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate my-python

/software/cellgen/team274/lr26/miniforge3/envs/my-python/bin/python /nfs/users/nfs_l/lr26/shells/sc_doublet.py
