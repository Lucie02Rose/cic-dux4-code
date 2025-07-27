#!/usr/bin/env bash
### This is the AlleleIntegrator wrapper script (used by the group) ###
CORES=52 # no of cpus
RAM="200G" # a lot of samples - a lot of memory
QUEUE="hugemem" # big memory queue
GROUP="team274" # group name

SCRIPT="Rscript /nfs/users/nfs_l/lr26/allele_check/genotypeCheck.R " # the wrapper runs this script
# this is a shared conda environment within the group (pre-Sanger datacentre incident)
conda activate alleleIntegrator
# LSF uses bsub to submit jobs
bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch126/casm/team274sb/lr26/outputs/allele-%J.output" -e "/lustre/scratch126/casm/team274sb/lr26/errors/allele-%J.error" \
"${SCRIPT}"
