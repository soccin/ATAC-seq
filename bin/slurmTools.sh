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
# bSync.sh runs in the background and is waited on. scancel and the time
# limit send SIGTERM to the batch shell only, and bash runs a trap only
# after its foreground child exits, so with bSync.sh in the foreground the
# caller's trap never ran before the SIGKILL KillWait (30 s) later. The
# wait builtin returns as soon as a trapped signal arrives. BSYNC_PID lets
# the trap stop the orphaned bSync.sh.
#
bSync() {
    bSync.sh "^$1\$" &
    BSYNC_PID=$!
    wait $BSYNC_PID
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

# bCheck JOBNAME [NJOBS]
#
# Call immediately after bSync JOBNAME. Aborts the calling script if any
# job with that name ended in a state other than COMPLETED. NJOBS, when
# given, is the number of jobs submitted under that name.
#
# Neither bSync nor sacct is taken at its word. bSync.sh reads a failed
# squeue call (controller timeout) as "no jobs left" and returns early,
# so a job sacct still shows as active sends bCheck back into bSync
# rather than counting as a failure. A failed sacct call, no rows, or
# fewer rows than NJOBS (accounting lags the queue) is retried every 30 s
# for up to ATAC_SACCT_WAIT seconds (default 1800) before it is fatal.
# Only a job in a final state other than COMPLETED fails the stage.
#
# Sleeps run in the background and are waited on, as in bSync, so a
# SIGTERM reaches the caller's trap at once.
#
bCheck() {
    local jobname=$1
    local njobs=${2:-1}
    local limit=${ATAC_SACCT_WAIT:-1800}
    local waited=0
    local rows nrows n_active n_bad

    while true; do
        if ! rows=$(sacctByName "$jobname"); then
            rows=""
            echo "bCheck: sacct failed for [$jobname]"
        fi
        nrows=$(echo "$rows" | awk 'NF > 0' | wc -l)

        if [ "$nrows" -lt "$njobs" ]; then
            if [ "$waited" -ge "$limit" ]; then
                echo
                echo "    FATAL: $nrows of $njobs accounting records for [$jobname]" \
                     "since $ATAC_START after ${waited}s"
                echo
                exit 1
            fi
            echo "bCheck: $nrows of $njobs accounting records for [$jobname]; retry in 30s"
            sleep 30 &
            wait $!
            waited=$((waited + 30))
            continue
        fi

        n_active=$(echo "$rows" \
            | awk -F'|' 'NF > 0 && $3 ~ /^(PENDING|RUNNING|COMPLETING|SUSPENDED|REQUEUED|CONFIGURING|RESIZING|SIGNALING|STAGE_OUT)/' \
            | wc -l)
        if [ "$n_active" -eq 0 ]; then
            break
        fi

        echo "bCheck: $n_active job(s) of [$jobname] still active in sacct; waiting again"
        bSync "$jobname"
        sleep 30 &
        wait $!
    done

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
