#!/bin/bash

#
# motifEnrich.sh -- motif enrichment of the differential peaks, as jobs
#
# usage, from the analysis directory, after R/diffAnalysisPairwise.R:
#   sbatch /path/to/ATAC-seq/motifEnrich.sh [-s] [-w WIDTH] [-m MINPEAKS] [RUNTAG]
#
# The comparisons are independent, so by default each is its own job.
# This job counts them (R/motifEnrich.R --count, which also checks the
# inputs), then submits an array job of one task per comparison
# (--comparison N, results saved in a directory on scratch) and a
# combine job (--combine) that runs when every task has COMPLETED and
# writes the output files; then it exits. If a task fails, Slurm
# cancels the combine job. The tasks and the combine job run
# bin/motifEnrichJob.sh with sbatch --wrap; their sizes are TASK_SBATCH
# and COMBINE_SBATCH below.
#
# With -s this job runs every comparison in turn, the only way before
# the per-comparison jobs were added, with the directives below: 8 cores, 128G and 6 h on cmobic_cpu with
# --qos=priority, as pipe.sh's long jobs; flags given to sbatch on the
# command line override them. The motif scan uses every core of the
# job (SLURM_CPUS_PER_TASK is the default of R/motifEnrich.R -c).
# Without -s this job only counts and submits, for about a minute.
#
# sbatch here is the ~/bin wrapper, as for pipe.sh. Slurm runs a copy
# of the batch script from its spool directory, so $0 cannot locate
# this checkout; the wrapper exports SBATCH_SCRIPT_DIR and SDIR is
# taken from that. Run directly (no sbatch), $0 is used.
#
# Outside a job, or with -h, it prints the sbatch line and how to
# choose the options and the job size (usage below). README.md, "Motif
# enrichment", has the same notes; keep the two the same.
# R/motifEnrich.R --help describes the method and the output.
#
#SBATCH -p cmobic_cpu
#SBATCH --qos=priority
#SBATCH -t 6:00:00
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -J MOTIF.ATAC
#SBATCH -o SLURM.MOTIF/%j.out

set -euo pipefail

SDIR=${SBATCH_SCRIPT_DIR:-$( cd "$( dirname "$0" )" && pwd -P )}

if [ ! -e "$SDIR/R/motifEnrich.R" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo "    submit with ~/bin/sbatch, which exports SBATCH_SCRIPT_DIR"
    echo
    exit 1
fi

# The size of each comparison's job and of the combine job. After the
# motif scan, monaLisa runs one worker per set (up, down and
# background), at most 3; 24 cores were no faster than 8, so 4
TASK_SBATCH=(-p cmobic_cpu --qos=priority -t 6:00:00 -c 4 --mem=96G)
# shellcheck disable=SC2054 # the comma separates the two partitions
COMBINE_SBATCH=(-p cpushort,cmobic_short -t 1:00:00 -c 1 --mem=32G)

usage() {
    cat <<EOF

usage: motifEnrich.sh [-s] [-w WIDTH] [-m MINPEAKS] [RUNTAG]

Submit from the analysis directory, after diffAnalysisPairwise.R,
with the same RUNTAG (default: none):

  sbatch $SDIR/motifEnrich.sh [-s] [-w WIDTH] [-m MINPEAKS] [RUNTAG]

Writes <PROJ>_<RUNTAG>_MotifEnrichment.xlsx, .pdf and .rds to
atacSeq/diff. By default each comparison is a job of its own and a
combine job writes the files; the results are the same as with -s.
The logs are in SLURM.MOTIF: <jobid>.out of this job, which prints
the ids of the others, <arrayid>_<N>.out of comparison N, and
<jobid>.out of the combine job. The files are complete when the
combine job is COMPLETED in sacct; if a comparison fails, Slurm
cancels the combine job.

-s
    Run every comparison in this job, one after another, the only
    way before the per-comparison jobs: one job and one log, but the
    time is the sum of the comparisons'.

-w WIDTH (default 501; odd, at least 51)
    The bp scanned for motifs around each peak's strongest summit, for
    the up, down and background peaks alike. 501 (summit +/- 250) is
    ArchR's peak width and in chromVAR's 250-500 range for ATAC, which
    it calls "fairly robust to the exact choice". HOMER's default,
    from TF ChIP-seq, is 200. Keep 501 unless you have a reason: 201
    (summit +/- 100) keeps closer to the summit; a wider window adds
    more flanking sequence. The scan time grows with the width.

-m MINPEAKS (default 50)
    An up or down set with fewer peaks is not tested; the Summary sheet
    says so. 50 is a judgment, not a standard: each set is tested for
    879 motifs, with BH. In Proj_18143 ED1 vs EV (33 up, 17 down) was
    not tested; of the sets tested, 244 peaks gave no motif with
    padj < 0.05 and 441 gave 92. Lower it to try a weak comparison and
    read its hits with care; raising it only skips more sets.

Time
    The time grows with the sequences scanned, about 220 per second
    with 8 cores (195 to 290 per comparison, Proj_18143). A comparison
    with at least one set tested scans all its tested peaks (the
    background too), so it takes about (tested peaks) / 220 s; one
    with no set tested is skipped. A comparison's tested peaks are the
    "N peaks" on its MA and volcano plots in
    <PROJ>_<RUNTAG>_DiffPeaks_V5.pdf. Proj_18143 scanned 6
    comparisons, 64,000 to 240,000 peaks each, 1.04 million in all.
    By default the run takes about as long as its largest comparison,
    with -s as long as the sum of them.

Job size, by default
    Each comparison: ${TASK_SBATCH[*]}
    Combine:         ${COMBINE_SBATCH[*]}
    These are TASK_SBATCH and COMBINE_SBATCH in motifEnrich.sh; flags
    given to sbatch apply only to this job, which needs none.
    Proj_18143_B (4 comparisons, about 230,000 tested peaks each):
    16 to 18 min per comparison on 4 cores (MaxRSS 54-56G with a 64G
    cap, so 96G for room), combine
    31 s (5G), 20 min in all; with -s 52 min (13.5 min per comparison
    on 8 cores, MaxRSS 66G). The results were identical.

Job size, with -s
    The script asks for -c 8 --mem=128G -t 6:00:00 (cmobic_cpu,
    --qos=priority). To change them, give the flags to sbatch before
    the script, e.g.
      sbatch -t 12:00:00 $SDIR/motifEnrich.sh -s
    -c      more cores do not shorten the run: Proj_18143 took 80 min
            with -c 8 and 78 min with -c 24 (fewer than 8 has not been
            measured).
    --mem   Proj_18143 reached 64G with -c 8 and 121G with -c 24 (it
            grows with the cores); 128G leaves room with -c 8.
    -t      the sum of the comparisons' times, with room; 6:00:00
            covers a run like Proj_18143.

The same notes: $SDIR/README.md,
"Motif enrichment", "Choosing the options and the job size".
The method and output: Rscript $SDIR/R/motifEnrich.R --help

EOF
}

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

SERIAL=no
R_ARGS=()
for arg in "$@"; do
    case "$arg" in
        -h|--help)
            usage
            exit 0
            ;;
        -s|--serial)
            SERIAL=yes
            ;;
        --count|--comparison*|--combine|--parts*)
            fatal "$arg is set by motifEnrich.sh itself"
            ;;
        *)
            R_ARGS+=("$arg")
            ;;
    esac
done

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo
    echo "    motifEnrich.sh must run as a job; submit it with sbatch"
    usage
    exit 1
fi
CTRL_JOBID=$SLURM_JOB_ID

. "$SDIR/bin/loadTools.sh"

echo "JOB=$SLURM_JOB_ID HOST=$(hostname) CPUS=${SLURM_CPUS_PER_TASK:-1} SERIAL=$SERIAL"
echo "SDIR=$SDIR"
echo "TMPDIR=$TMPDIR"
echo "PWD=$(pwd -P)"
date

if [[ "$SERIAL" == "yes" ]]; then
    Rscript "$SDIR/R/motifEnrich.R" "${R_ARGS[@]}"
    date
    exit 0
fi

#
# One job per comparison. --count checks the inputs and the options
# before anything is submitted.
#

if ! COUNT_OUT=$(Rscript "$SDIR/R/motifEnrich.R" --count "${R_ARGS[@]}"); then
    echo "$COUNT_OUT"
    fatal "R/motifEnrich.R --count failed"
fi
echo "$COUNT_OUT"
NCOMP=$(sed -n 's/^COMPARISONS=//p' <<< "$COUNT_OUT")
if ! [[ "$NCOMP" =~ ^[1-9][0-9]*$ ]]; then
    fatal "no number of comparisons from R/motifEnrich.R --count [$NCOMP]"
fi

# On the shared scratch, which every node sees; one directory per run.
# The combine job removes it; it is kept when a job fails.
PARTS=${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}/motif.$CTRL_JOBID
if ! mkdir -p "$PARTS"; then
    fatal "cannot create $PARTS"
fi

# Nothing of this job's Slurm environment may reach the jobs it
# submits, as in pipe.sh: SLURM_* but SLURM_CONF, and every SBATCH_*
# (SBATCH_SCRIPT_DIR too; SDIR is set)
for V in ${!SLURM_@}; do
    if [ "$V" != "SLURM_CONF" ]; then
        unset "$V"
    fi
done
for V in ${!SBATCH_@}; do
    unset "$V"
done

# --wrap runs the command with sh; SLURM_ARRAY_TASK_ID is expanded in
# the task, everything else is quoted here
JOB_CMD=$(printf '%q ' bash "$SDIR/bin/motifEnrichJob.sh")
COMMON_ARGS=$(printf '%q ' --parts "$PARTS" "${R_ARGS[@]}")
TASK_WRAP="$JOB_CMD --comparison \"\$SLURM_ARRAY_TASK_ID\" $COMMON_ARGS"
COMBINE_WRAP="$JOB_CMD --combine $COMMON_ARGS"

if ! ARRAY_OUT=$(sbatch --parsable "${TASK_SBATCH[@]}" --array=1-"$NCOMP" \
        -J MOTIF.CMP."$CTRL_JOBID" -o "SLURM.MOTIF/%A_%a.out" \
        --wrap "$TASK_WRAP"); then
    fatal "cannot submit the comparison jobs"
fi
ARRAY_JOBID=${ARRAY_OUT%%;*}

if ! COMBINE_OUT=$(sbatch --parsable "${COMBINE_SBATCH[@]}" \
        --dependency=afterok:"$ARRAY_JOBID" --kill-on-invalid-dep=yes \
        -J MOTIF.COMBINE."$CTRL_JOBID" -o "SLURM.MOTIF/%j.out" \
        --wrap "$COMBINE_WRAP"); then
    if ! scancel "$ARRAY_JOBID"; then
        echo "    cannot cancel the comparison jobs $ARRAY_JOBID; cancel them by hand"
    fi
    fatal "cannot submit the combine job; comparison jobs $ARRAY_JOBID cancelled"
fi
COMBINE_JOBID=${COMBINE_OUT%%;*}

echo
echo "Comparisons: $NCOMP, array job $ARRAY_JOBID (${TASK_SBATCH[*]})"
echo "    logs SLURM.MOTIF/${ARRAY_JOBID}_<N>.out, results in $PARTS"
echo "Combine: job $COMBINE_JOBID (${COMBINE_SBATCH[*]}), after every comparison"
echo "    log SLURM.MOTIF/$COMBINE_JOBID.out"
echo "The files are complete when job $COMBINE_JOBID is COMPLETED:"
echo "    sacct -j $ARRAY_JOBID,$COMBINE_JOBID -o JobID,JobName%24,State,Elapsed,MaxRSS"
date
