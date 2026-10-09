#!/bin/bash
#
# slurmTools.sh -- helpers for the control-job style pipe.sh on IRIS/Slurm
#
# Usage: source slurmTools.sh
#
# Jobs are submitted with ~/bin/bsub (the SchedMD openlava shim with local
# patches: -W in sbatch time format, -M for memory, -o DIR/ for per-job
# logs) and waited on by job name with ~/bin/bSync.sh.
#
# Slurm writes nothing of its own into a job's log, so a job that was
# OOM-killed, timed out or exited non-zero leaves a log that can look like
# a good one. "Did it work" is therefore answered by sacct, matched on
# State: anything but COMPLETED is a failure. Match State, never ExitCode:
# an OOM kill reports State=OUT_OF_MEMORY with ExitCode=0:125.
#
# Job names embed the control job's $$ and pids recycle, so every sacct
# query is bounded by ATAC_START. pipe.sh exports it before the first
# submission; standalone use falls back to the last hour.
#

ATAC_START=${ATAC_START:-$(date -d '-1 hour' +%Y-%m-%dT%H:%M:%S)}

# bSync JOBNAME
#
# Block until no job of this user with exactly that name is pending or
# running. bSync.sh takes a regex, so anchor it: unanchored, JOBNAME would
# also match JOBNAME<more digits> from a concurrent run.
#
bSync() {
    bSync.sh "^$1\$"
}

# sacctByName JOBNAME
#
# One line per job with that name since ATAC_START:
#   JobID|JobName|State|ExitCode|Elapsed
#
sacctByName() {
    sacct -X -n -P -u "$USER" -S "$ATAC_START" --name="$1" \
        --format=JobID,JobName,State,ExitCode,Elapsed
}

# bCheck JOBNAME
#
# Call immediately after bSync JOBNAME. Aborts the calling script if any
# job with that name did not reach COMPLETED. Accounting can lag the
# queue by a few seconds, so an empty answer, or one that still shows a
# job in flight, is retried before it is treated as fatal.
#
bCheck() {
    local jobname=$1
    local rows n_active n_bad try

    for try in 1 2 3 4 5 6; do
        rows=$(sacctByName "$jobname")
        n_active=$(echo "$rows" \
            | awk -F'|' 'NF > 0 && $3 ~ /^(PENDING|RUNNING|COMPLETING|SUSPENDED|REQUEUED|CONFIGURING)/' \
            | wc -l)
        if [ -n "$rows" ] && [ "$n_active" == "0" ]; then
            break
        fi
        echo "bCheck: waiting for accounting records for [$jobname] (try $try)"
        sleep 10
    done

    if [ -z "$rows" ]; then
        echo
        echo "    FATAL: no accounting records for [$jobname] since $ATAC_START"
        echo
        exit 1
    fi

    n_bad=$(echo "$rows" | awk -F'|' 'NF > 0 && $3 !~ /^COMPLETED/' | wc -l)
    if [ "$n_bad" -gt 0 ]; then
        echo
        echo "    FATAL: $n_bad job(s) failed for [$jobname]"
        echo "$rows" \
            | awk -F'|' '$3 !~ /^COMPLETED/ {printf "    %-10s %-14s %-8s %-10s %s\n", $1, $3, $4, $5, $2}'
        echo
        exit 1
    fi

    echo "bCheck: [$jobname] $(echo "$rows" | wc -l) job(s) COMPLETED"
}

# bCheckAll REGEX
#
# Print every job since ATAC_START whose name matches REGEX and that is
# not COMPLETED. Prints nothing when all is well. This is the Slurm
# version of the old "find LSF* | parseLSF.py | fgrep -v Successfully"
# sweep at the end of pipe.sh.
#
bCheckAll() {
    sacct -X -n -P -u "$USER" -S "$ATAC_START" \
        --format=JobID,JobName,State,ExitCode,Elapsed \
        | awk -F'|' -v re="$1" \
            '$2 ~ re && $3 !~ /^COMPLETED/ {printf "    %-10s %-14s %-8s %-10s %s\n", $1, $3, $4, $5, $2}'
}

# bKill REGEX
#
# scancel every queued or running job of this user whose name matches
# REGEX. pipe.sh calls it from its EXIT trap so an aborted run does not
# leave the rest of the stage in the queue.
#
bKill() {
    squeue -h -u "$USER" -o "%i %j" \
        | awk -v re="$1" '$2 ~ re {print $1}' \
        | xargs -r scancel
}
