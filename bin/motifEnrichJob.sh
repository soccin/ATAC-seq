#!/bin/bash

#
# motifEnrichJob.sh -- run R/motifEnrich.R in a job with the tools of
# bin/loadTools.sh
#
# usage: bash $SDIR/bin/motifEnrichJob.sh MOTIFENRICH_R_ARGS ...
#
# motifEnrich.sh (without -s) submits it with sbatch --wrap, once per
# comparison (--comparison N --parts DIR) and once to combine them
# (--combine --parts DIR). It is run by its path, not as a batch
# script, so $0 locates this checkout. It exits with R's status, so a
# failed comparison is FAILED in sacct and cancels the combine job.
#

set -euo pipefail

SDIR=$( cd "$( dirname "$0" )/.." && pwd -P )

if [ ! -e "$SDIR/R/motifEnrich.R" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo
    exit 1
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo
    echo "    motifEnrichJob.sh must run as a job; motifEnrich.sh submits it"
    echo
    exit 1
fi

. "$SDIR/bin/loadTools.sh"

echo "JOB=$SLURM_JOB_ID ARRAY=${SLURM_ARRAY_JOB_ID:-none}_${SLURM_ARRAY_TASK_ID:-none}" \
    "HOST=$(hostname) CPUS=${SLURM_CPUS_PER_TASK:-1}"
echo "SDIR=$SDIR"
echo "TMPDIR=$TMPDIR"
echo "PWD=$(pwd -P)"
echo "ARGS=$*"
date
Rscript "$SDIR/R/motifEnrich.R" "$@"
date
