#!/bin/bash

#
# motifEnrichJob.sh -- run R/motifEnrich.R in a job with the tools of
# bin/loadTools.sh
#
# usage: bash $SDIR/bin/motifEnrichJob.sh [--status RECORD] MOTIFENRICH_R_ARGS ...
#
# motifEnrich.sh (without -s) submits it with sbatch --wrap, once per
# comparison (--comparison N --parts DIR) and once to combine them
# (--combine --parts DIR). It is run by its path, not as a batch
# script, so $0 locates this checkout. It exits with R's status, so a
# failed comparison is FAILED in sacct and cancels the combine job.
# With --status, it adds STATUS=FAILED to the run's record
# (bin/motifStatus.sh) when it fails, and STATUS=COMPLETED when it is
# the combine job and has written the files.
#

set -euo pipefail

SDIR=$( cd "$( dirname "$0" )/.." && pwd -P )

if [ ! -e "$SDIR/R/motifEnrich.R" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo
    exit 1
fi

STATUSF=""
if [[ "${1:-}" == "--status" ]]; then
    STATUSF=$2
    shift 2
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo
    echo "    motifEnrichJob.sh must run as a job; motifEnrich.sh submits it"
    echo
    exit 1
fi

COMBINE=no
for arg in "$@"; do
    case "$arg" in
        --combine)
            COMBINE=yes
            ;;
    esac
done
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    JOB_LOG=SLURM.MOTIF/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
else
    JOB_LOG=SLURM.MOTIF/$SLURM_JOB_ID.out
fi

recordExit() {
    local rc=$?
    if [ "$rc" != "0" ]; then
        {
            echo "STATUS=FAILED"
            echo "MESSAGE=job $SLURM_JOB_ID exited $rc, log $JOB_LOG"
            echo "FINISHED=$(date '+%Y-%m-%d %H:%M:%S')"
        } >> "$STATUSF"
    elif [[ "$COMBINE" == "yes" ]]; then
        {
            echo "STATUS=COMPLETED"
            echo "FINISHED=$(date '+%Y-%m-%d %H:%M:%S')"
        } >> "$STATUSF"
    fi
    exit "$rc"
}
if [[ -n "$STATUSF" ]]; then
    trap recordExit EXIT
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
