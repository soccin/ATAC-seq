#!/bin/bash

#
# checkRun.sh -- did the ATAC-seq run in an analysis directory work?
#
# usage:
#   checkRun.sh [ANALYSIS_DIR]      default: the current directory
#
# exit status:
#   0   the run completed and every deliverable is present
#   1   the run failed
#   2   the run is still going
#   3   unknown: sacct failed, so the state of the run cannot be read now
#
# Slurm writes nothing of its own into a job's -o file, so no log can
# answer this by itself. pipe.sh leaves three records and this script
# reads them together with sacct:
#
#   00.RUNSTATUS.txt     STATUS=RUNNING|COMPLETED|FAILED, the STAGE it
#                        reached and the control job id. COMPLETED is
#                        written only after the per-sample deliverables
#                        check passes.
#
#   SLURM.CTRL/jobs.tsv  every stage job: job id, stage, sample, log.
#
#   #ATAC_EXIT=<rc>      the trailer bin/runStage.sh appends to every
#                        stage log. Used only when sacct no longer has
#                        the job.
#
# STATUS=RUNNING is not taken at its word: a control job killed outright
# (memory cap, node failure) cannot update the file, so the control job's
# own state in sacct decides whether the run is really still going.
#

usage () {
    echo
    echo "usage: checkRun.sh [ANALYSIS_DIR]"
    echo
    echo "  exit: 0 worked, 1 failed, 2 still running, 3 unknown (sacct failed)"
    echo
    exit 1
}

case "$1" in
    -h|--help)
        usage
        ;;
esac

DIR=${1:-.}
STATUSF=$DIR/00.RUNSTATUS.txt
JOBS=$DIR/SLURM.CTRL/jobs.tsv

if [ ! -e "$STATUSF" ]; then
    echo "ATAC STATUS: NONE   no $STATUSF; pipe.sh has not started here"
    exit 1
fi

getKey () {
    sed -n "s/^$1=//p" "$STATUSF" | tail -1
}

STATUS=$(getKey STATUS)
STAGE=$(getKey STAGE)
RC=$(getKey RC)
MESSAGE=$(getKey MESSAGE)
STARTED=$(getKey STARTED)
FINISHED=$(getKey FINISHED)
VERSION=$(getKey VERSION)
CTRL_JOBID=$(getKey CTRL_JOBID)
CTRL_LOG=$(getKey CTRL_LOG)
HOST=$(getKey HOST)
PID=$(getKey PID)
SAMPLES=$(getKey SAMPLES)

ACTIVE_RE='^(PENDING|RUNNING|SUSPENDED|COMPLETING|CONFIGURING|REQUEUED|RESIZING|SIGNALING|STAGE_OUT)$'

##
# Stage jobs: sacct is the authority; the log trailer is the fallback for
# jobs sacct has purged. "CANCELLED by 12345" keeps only the first word.

NJOBS=0
JOBREPORT=""
SACCT_FAILED=0
if [ -e "$JOBS" ]; then
    IDS=$(awk -F'\t' '/^#/ {next} NF > 0 {print $1}' "$JOBS" | tr '\n' ',' | sed 's/,$//')
    if [ -n "$IDS" ]; then
        if SACCT=$(sacct -X -n -P -j "$IDS" --format=JobID,State,ExitCode,Elapsed 2>/dev/null); then
            SACCT=$(printf '%s\n' "$SACCT" | tr '|' '\t')
        else
            SACCT=""
            SACCT_FAILED=1
        fi
        JOBREPORT=$(awk -F'\t' -v dir="$DIR" -v active="$ACTIVE_RE" '

            function trailer(logf,   cmd, line, rc) {
                rc = ""
                cmd = "tail -20 \"" logf "\" 2>/dev/null"
                while ((cmd | getline line) > 0) {
                    if (match(line, /#ATAC_EXIT=[0-9]+/)) {
                        rc = substr(line, RSTART + 11, RLENGTH - 11)
                    }
                }
                close(cmd)
                return rc
            }

            NR == FNR {
                split($2, s, " ")
                state[$1] = s[1]
                code[$1] = $3
                secs[$1] = $4
                next
            }

            /^#/ { next }
            NF == 0 { next }

            {
                id = $1
                stage = $2
                sample = $3
                logf = dir "/" $4

                st = state[id]
                if (st == "") {
                    rc = trailer(logf)
                    if (rc == "0") {
                        st = "COMPLETED"
                    } else if (rc != "") {
                        st = "FAILED"
                        code[id] = rc ":0"
                    } else {
                        st = "NO-RECORD"
                    }
                    secs[id] = "-"
                }

                njob++
                if (st == "COMPLETED") {
                    nok++
                    next
                }
                if (st ~ active) {
                    nactive++
                    act[nactive] = sprintf("   %-10s %-12s %-10s %s", id, st, stage, sample)
                    next
                }
                #
                # CANCELLED is almost always bKill clearing the rest of
                # the run after a real failure, so list it apart.
                #
                if (st == "CANCELLED") {
                    ncancel++
                    canc[ncancel] = sprintf("   %-10s %-12s %-10s %s", id, st, stage, sample)
                    next
                }
                nfail++
                bad[nfail] = sprintf("   %-10s %-14s %-8s %-10s %-10s %s\n              %s", \
                                     id, st, code[id], secs[id], stage, sample, logf)
            }

            END {
                printf("JOBS\t%d\t%d\t%d\t%d\t%d\n", njob, nok, nfail, ncancel, nactive)
                if (nfail) {
                    print "   ---- failed"
                    for (i = 1; i <= nfail; i++) { print bad[i] }
                }
                if (ncancel) {
                    print "   ---- cancelled"
                    for (i = 1; i <= ncancel; i++) { print canc[i] }
                }
                if (nactive) {
                    print "   ---- queued or running"
                    for (i = 1; i <= nactive; i++) { print act[i] }
                }
            }

        ' <(printf '%s\n' "$SACCT") "$JOBS")
    fi
fi

COUNTS=$(echo "$JOBREPORT" | awk -F'\t' '$1 == "JOBS"')
NJOBS=$(echo "$COUNTS" | cut -f2)
NOK=$(echo "$COUNTS" | cut -f3)
NFAIL=$(echo "$COUNTS" | cut -f4)
NCANCEL=$(echo "$COUNTS" | cut -f5)
NACTIVE=$(echo "$COUNTS" | cut -f6)
NJOBS=${NJOBS:-0}
NOK=${NOK:-0}
NFAIL=${NFAIL:-0}
NCANCEL=${NCANCEL:-0}
NACTIVE=${NACTIVE:-0}
JOBLIST=$(echo "$JOBREPORT" | awk -F'\t' '$1 != "JOBS"')

##
# Control job

CTRL_STATE=""
if [ "$CTRL_JOBID" != "none" ] && [ -n "$CTRL_JOBID" ]; then
    if CTRL_SACCT=$(sacct -X -n -P -j "$CTRL_JOBID" --format=State 2>/dev/null); then
        CTRL_STATE=$(printf '%s\n' "$CTRL_SACCT" | head -1 | awk '{print $1}')
    else
        SACCT_FAILED=1
    fi
fi

NOTE=""
case "$STATUS" in

    COMPLETED)
        if [ "$NFAIL" == "0" ] && [ "$NCANCEL" == "0" ] && [ "$NACTIVE" == "0" ]; then
            VERDICT=OK
            EXIT=0
        else
            VERDICT=FAILED
            EXIT=1
            NOTE="status file says COMPLETED but not every stage job did"
        fi
        ;;

    FAILED)
        VERDICT=FAILED
        EXIT=1
        ;;

    RUNNING)
        if [ "$CTRL_JOBID" == "none" ]; then
            #
            # pipe.sh run by hand, not under Slurm: only the pid can say.
            #
            if [ "$HOST" != "$(hostname)" ]; then
                VERDICT=RUNNING
                EXIT=2
                NOTE="pipe.sh ran outside Slurm on $HOST; run checkRun.sh there to confirm"
            elif kill -0 "$PID" 2>/dev/null; then
                VERDICT=RUNNING
                EXIT=2
            else
                VERDICT=FAILED
                EXIT=1
                NOTE="pipe.sh (pid $PID) is gone but never set a final status"
            fi
        elif [[ "$CTRL_STATE" =~ $ACTIVE_RE ]]; then
            VERDICT=RUNNING
            EXIT=2
        else
            VERDICT=FAILED
            EXIT=1
            NOTE="control job $CTRL_JOBID ended (${CTRL_STATE:-no sacct record}) without setting a final status"
        fi
        ;;

    *)
        VERDICT=FAILED
        EXIT=1
        NOTE="unrecognized STATUS [$STATUS] in $STATUSF"
        ;;

esac

#
# A failed sacct call (controller timeout) leaves stage jobs with no state
# and the control job with none either, which reads as FAILED. Only a
# STATUS=FAILED written by pipe.sh, or an OK verdict, stands without sacct.
#
if [ "$SACCT_FAILED" == "1" ] && [ "$EXIT" != "0" ] \
       && { [ "$STATUS" == "RUNNING" ] || [ "$STATUS" == "COMPLETED" ]; }; then
    VERDICT=UNKNOWN
    EXIT=3
    NOTE="sacct failed; cannot tell whether the run is going, worked or failed. Run checkRun.sh again."
    JOBLIST=""
fi

echo "ATAC STATUS: $VERDICT   $(date '+%Y-%m-%d %H:%M:%S')"
echo "   dir        $(cd "$DIR" && pwd)"
echo "   version    $VERSION"
echo "   control    $CTRL_JOBID ${CTRL_STATE:+($CTRL_STATE)}  $CTRL_LOG"
echo "   started    $STARTED"
if [ -n "$FINISHED" ]; then
    echo "   finished   $FINISHED"
fi
echo "   stage      $STAGE${RC:+  rc=$RC}"
echo "   samples    $SAMPLES"
if [ "$VERDICT" == "UNKNOWN" ]; then
    echo "   jobs       $NJOBS: states not available (sacct failed)"
else
    echo "   jobs       $NJOBS: $NOK ok, $NFAIL failed, $NCANCEL cancelled, $NACTIVE queued/running"
fi
if [ -n "$MESSAGE" ]; then
    echo "   message    $MESSAGE"
fi
if [ -n "$NOTE" ]; then
    echo "   NOTE       $NOTE"
fi
if [ -n "$JOBLIST" ]; then
    echo "$JOBLIST"
fi

exit $EXIT
