#!/bin/bash

#
# motifStatus.sh -- the state of a motif enrichment run, for motifEnrich.sh
#
# Sourced by motifEnrich.sh. motifStatus STATUSF OUTBASE reads the
# run's record STATUSF (SLURM.MOTIF/00.MOTIFSTATUS[_RUNTAG].txt) and
# the jobs it names in sacct, prints a report and sets
#
#   MOTIF_VERDICT  OK      the output files are complete
#                  FAILED  a job failed, or the run ended without them
#                  RUNNING a job of the run is queued or running
#                  UNKNOWN sacct failed, or the files exist with no record
#                  NONE    no record and no output files
#   MOTIF_EXIT     0, 1, 2, 3 and 4 for these
#   MOTIF_ACTIVE   the number of the run's jobs queued or running
#   MOTIF_HOW      SACCT when sacct could not be read, NORECORD when
#                  the output files exist but there is no record
#
# The record holds KEY=VALUE lines, appended; the last of a key counts.
# motifEnrich.sh writes STATUS=RUNNING when it starts, SUBMITTED once
# the comparison and combine jobs are queued (COMPLETED at the end with
# -s), FAILED from its EXIT trap; bin/motifEnrichJob.sh adds COMPLETED
# when the combine job has written the files and FAILED when a job
# fails. A job killed outright (memory cap, time limit) writes nothing,
# so sacct decides whenever the record does not say COMPLETED or FAILED.
#

# The MOTIF_* results are read by motifEnrich.sh, which sources this
# shellcheck disable=SC2034
MOTIF_ACTIVE_RE='^(PENDING|RUNNING|SUSPENDED|COMPLETING|CONFIGURING|REQUEUED|RESIZING|SIGNALING|STAGE_OUT)$'

motifKey() {
    sed -n "s/^$2=//p" "$1" | tail -1
}

motifStatus() {
    local statusf=$1
    local outbase=$2
    local ext missing nmissing status mode ctrl ctrl_log array combine
    local ids sacct_out sacct_failed summary sumline joblist note
    local ctrl_state combine_state tasks nactive nbad args

    MOTIF_VERDICT=""
    MOTIF_EXIT=""
    MOTIF_ACTIVE=0
    MOTIF_HOW=""
    note=""
    joblist=""

    missing=""
    nmissing=0
    for ext in xlsx pdf rds; do
        if [ ! -s "$outbase.$ext" ]; then
            missing="$missing $outbase.$ext"
            nmissing=$((nmissing + 1))
        fi
    done

    if [ ! -e "$statusf" ]; then
        if [ "$nmissing" -lt 3 ]; then
            MOTIF_VERDICT=UNKNOWN
            MOTIF_EXIT=3
            MOTIF_HOW=NORECORD
            note="output files exist but there is no record: made before runs were recorded, or by hand"
        else
            MOTIF_VERDICT=NONE
            MOTIF_EXIT=4
        fi
        echo "MOTIF STATUS: $MOTIF_VERDICT   $(date '+%Y-%m-%d %H:%M:%S')"
        echo "   dir        $(pwd -P)"
        echo "   record     $statusf (none)"
        echo "   outputs    $outbase.{xlsx,pdf,rds}"
        if [ "$nmissing" -gt 0 ] && [ "$nmissing" -lt 3 ]; then
            echo "   missing   $missing"
        fi
        if [ -n "$note" ]; then
            echo "   NOTE       $note"
        fi
        return 0
    fi

    status=$(motifKey "$statusf" STATUS)
    mode=$(motifKey "$statusf" MODE)
    ctrl=$(motifKey "$statusf" CTRL_JOBID)
    ctrl_log=$(motifKey "$statusf" CTRL_LOG)
    array=$(motifKey "$statusf" ARRAY_JOBID)
    combine=$(motifKey "$statusf" COMBINE_JOBID)
    if [ "$ctrl" == "none" ]; then
        ctrl=""
    fi

    ids=$(echo "$ctrl,$array,$combine" | sed 's/,,*/,/g; s/^,//; s/,$//')
    sacct_failed=0
    sacct_out=""
    if [ -n "$ids" ]; then
        if ! sacct_out=$(sacct -X -n -P -j "$ids" --format=JobID,JobName,State,ExitCode,Elapsed 2>/dev/null); then
            sacct_failed=1
            sacct_out=""
        fi
    fi

    #
    # One line of states, then a line per job neither COMPLETED nor
    # active. "CANCELLED by 12345" keeps only the first word.
    #
    summary=$(awk -F'|' -v ctrl="$ctrl" -v arr="$array" -v comb="$combine" \
                  -v ctrllog="$ctrl_log" -v active="$MOTIF_ACTIVE_RE" '

        NF == 0 { next }

        {
            id = $1
            split($3, s, " ")
            st = s[1]
            if (id == ctrl) {
                cst = st
                logf = ctrllog
            } else if (id == comb) {
                bst = st
                logf = "SLURM.MOTIF/" id ".out"
            } else if (arr != "" && index(id, arr "_") == 1) {
                ntask++
                tst[st]++
                logf = "SLURM.MOTIF/" id ".out"
            } else {
                next
            }
            if (st ~ active) {
                nact++
            } else if (st != "COMPLETED") {
                nbad++
                bad[nbad] = sprintf("   %-14s %-26s %-12s %-6s %-9s %s", \
                                    id, $2, st, $4, $5, logf)
            }
        }

        END {
            tasks = ""
            for (st in tst) {
                tasks = tasks (tasks == "" ? "" : ", ") tst[st] " " st
            }
            printf("SUM\t%s\t%s\t%s\t%d\t%d\n", cst, bst, tasks, nact, nbad)
            for (i = 1; i <= nbad; i++) {
                print bad[i]
            }
        }

    ' <<< "$sacct_out")

    sumline=$(echo "$summary" | head -1)
    joblist=$(echo "$summary" | tail -n +2)
    ctrl_state=$(echo "$sumline" | cut -f2)
    combine_state=$(echo "$sumline" | cut -f3)
    tasks=$(echo "$sumline" | cut -f4)
    nactive=$(echo "$sumline" | cut -f5)
    nbad=$(echo "$sumline" | cut -f6)
    MOTIF_ACTIVE=$nactive

    if [ "$sacct_failed" == "1" ]; then
        MOTIF_HOW=SACCT
        joblist=""
        if [ "$status" == "COMPLETED" ] && [ "$nmissing" == "0" ]; then
            MOTIF_VERDICT=OK
            MOTIF_EXIT=0
        elif [ "$status" == "FAILED" ]; then
            MOTIF_VERDICT=FAILED
            MOTIF_EXIT=1
            note="sacct failed; jobs of this run may still be queued or running"
        else
            MOTIF_VERDICT=UNKNOWN
            MOTIF_EXIT=3
            note="sacct failed; cannot tell whether the run is going, worked or failed. Check again."
        fi
    elif [ "$nactive" -gt 0 ]; then
        if [ "$status" == "FAILED" ] || [ "$nbad" -gt 0 ]; then
            MOTIF_VERDICT=FAILED
            MOTIF_EXIT=1
            note="a job failed, so the run cannot complete; $nactive job(s) still queued or running"
        else
            MOTIF_VERDICT=RUNNING
            MOTIF_EXIT=2
        fi
    elif [ "$status" == "COMPLETED" ]; then
        if [ "$nmissing" == "0" ]; then
            MOTIF_VERDICT=OK
            MOTIF_EXIT=0
        else
            MOTIF_VERDICT=FAILED
            MOTIF_EXIT=1
            note="the record says COMPLETED, but output files are missing"
        fi
    elif [ "$status" == "FAILED" ]; then
        MOTIF_VERDICT=FAILED
        MOTIF_EXIT=1
    elif [ "$combine_state" == "COMPLETED" ] && [ "$nmissing" == "0" ]; then
        MOTIF_VERDICT=OK
        MOTIF_EXIT=0
        note="the combine job COMPLETED but did not update the record"
    else
        MOTIF_VERDICT=FAILED
        MOTIF_EXIT=1
        note="no job of the run is queued or running, and it never recorded completion (STATUS=$status)"
    fi

    echo "MOTIF STATUS: $MOTIF_VERDICT   $(date '+%Y-%m-%d %H:%M:%S')"
    echo "   dir        $(pwd -P)"
    echo "   record     $statusf"
    args=$(motifKey "$statusf" ARGS)
    echo "   args       ${args:-(none)}"
    echo "   mode       $mode"
    echo "   started    $(motifKey "$statusf" STARTED)"
    if [ -n "$(motifKey "$statusf" FINISHED)" ]; then
        echo "   finished   $(motifKey "$statusf" FINISHED)"
    fi
    if [ -n "$(motifKey "$statusf" CTRL_JOBID)" ]; then
        echo "   control    $(motifKey "$statusf" CTRL_JOBID)${ctrl_state:+ ($ctrl_state)}  $ctrl_log"
    fi
    if [ -n "$array" ]; then
        echo "   tasks      $array, $(motifKey "$statusf" NCOMP) comparisons${tasks:+: $tasks}  SLURM.MOTIF/${array}_<N>.out"
    fi
    if [ -n "$combine" ]; then
        echo "   combine    $combine ${combine_state:+($combine_state)}  SLURM.MOTIF/$combine.out"
    fi
    echo "   outputs    $outbase.{xlsx,pdf,rds}"
    if [ "$nmissing" -gt 0 ] && [ "$MOTIF_VERDICT" != "RUNNING" ]; then
        echo "   missing   $missing"
    fi
    if [ -n "$(motifKey "$statusf" MESSAGE)" ]; then
        echo "   message    $(motifKey "$statusf" MESSAGE)"
    fi
    if [ -n "$note" ]; then
        echo "   NOTE       $note"
    fi
    if [ -n "$joblist" ]; then
        echo "   ---- not completed"
        echo "$joblist"
    fi
    return 0
}
