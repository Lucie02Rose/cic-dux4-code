#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 32
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q long
#BSUB -J dmp
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/outputs/%J-deepvar-mom.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/errors/%J-deepvar-mom.err
### note that this script was used pre-farmageddon (hence casm/team274sb)

### load teh farm singulatiry module ###
module load ISG/singularity/3.11.4

### set and export caches for the singulaty ###
export SINGULARITY_CACHEDIR=/lustre/scratch126/casm/team274sb/lr26/singularity
export SINGULARITY_TMPDIR=/lustre/scratch126/casm/team274sb/lr26/singularity/tmp

### define the paths used by the container ### 
SINGULARITY_IMG="/lustre/scratch126/casm/team274sb/lr26/singularity/pepper_margin_deepvariant_r0.8.sif"
BAM_FILE="/lustre/scratch126/casm/team274sb/lr26/PacBio-aligned/mom_1B02_hifi_reads.bam"
REFERENCE="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
OUTPUT_DIR="/lustre/scratch126/casm/team274sb/lr26/PacBio-deepvariant-mom"

### make the output directory ###
mkdir -p $OUTPUT_DIR

### run the pepper margin deepvariant container with singularity according to instructions on their github ###
singularity exec --bind /lustre/scratch126/casm/team274sb/lr26:/mnt \
    "${SINGULARITY_IMG}" \
    run_pepper_margin_deepvariant call_variant \
    -b "/mnt/PacBio-aligned/mom_1B02_hifi_reads.bam" \
    -f "/mnt/T2T/chm13v2.0.fa" \
    -o "/mnt/PacBio-deepvariant-mom" \
    -p "mom_output" \
    -t 32 \
    --hifi
