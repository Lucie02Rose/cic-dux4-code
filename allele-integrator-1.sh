#!/usr/bin/env bash

CORES=52
RAM="200G"
QUEUE="hugemem"
GROUP="team274"

SCRIPT="Rscript /nfs/users/nfs_l/lr26/allele_check/genotypeCheck.R "

#conda activate alleleIntegrator

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch126/casm/team274sb/lr26/output_logs/%J.output" -e "/lustre/scratch126/casm/team274sb/lr26/error_logs/%J.error" \
"${SCRIPT}"
