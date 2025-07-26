#!/usr/bin/env bash
# This is the AlleleIntegrator wrapper script (used by the group)
CORES=52
RAM="200G"
QUEUE="hugemem"
GROUP="team274"

SCRIPT="Rscript /nfs/users/nfs_l/lr26/allele_check/genotypeCheck.R "
# this is a shared conda environment within the group (pre-Sanger datacentre incident)
conda activate alleleIntegrator

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch126/casm/team274sb/lr26/outputs/allele-%J.output" -e "/lustre/scratch126/casm/team274sb/lr26/errors/allele-%J.error" \
"${SCRIPT}"
