#!/bin/bash

#
# Control job for IRIS/Slurm. Run it from the analysis directory; every
# stage is fanned out with ~/bin/bsub onto the short partitions and the
# control job blocks between stages with bSync/bCheck. The control job
# itself can run for hours on full-size BAMs, so the directives below put
# it on the long partition with --qos=priority. Flags given to sbatch on
# the command line override them.
#
# CMD:
#    mkdir -p SLURM.CTRL
#    sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]
#
# Then, at any time, from the same directory:
#
#    /path/to/ATAC-seq/bin/checkRun.sh
#
# which exits 0 if the run worked, 1 if it failed and 2 if it is still
# running.
#
# sbatch here is the ~/bin wrapper. Slurm runs a copy of the batch script
# from its spool directory, so $0 cannot locate this checkout; the wrapper
# exports SBATCH_SCRIPT_DIR and SDIR is taken from that. Running pipe.sh
# directly (no sbatch) falls back to $0.
#
#SBATCH -p cmobic_cpu
#SBATCH --qos=priority
#SBATCH -t 3-00:00:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -J CTRL.ATAC
#SBATCH -o SLURM.CTRL/%j.out

set -e

SDIR=${SBATCH_SCRIPT_DIR:-$( cd "$( dirname "$0" )" && pwd )}

if [ ! -e "$SDIR/postMapBamProcessing_ATACSeq.sh" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo "    submit with ~/bin/sbatch, which exports SBATCH_SCRIPT_DIR"
    echo
    exit 1
fi

source $SDIR/bin/loadTools.sh
source $SDIR/bin/slurmTools.sh

if [ ! -e "$SDIR/venv" ]; then
    echo
    echo "   Need to install macs2"
    echo "   Info in README"
    echo
    exit 1
fi

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
PIPENAME="ATAC-Seq"

##
# Process command args

MAPQ=10

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -q|--mapq)
      MAPQ="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

TAG=q$PIPENAME

COMMAND_LINE=$*

function usage {
    echo
    echo "usage: $PIPENAME/pipe.sh [-q MAPQ] BAM1 [BAM2 ... BAMN]"
    echo "version=$SCRIPT_VERSION"
    echo ""
    echo "Default MAPQ==$MAPQ"
    echo
    exit 1
}

if [ "$#" -lt "1" ]; then
    usage
fi

echo "MAPQ = ${MAPQ}"

BAMS=$*
echo SDIR=$SDIR
echo BAMS=$BAMS
echo VERSION=$SCRIPT_VERSION

#
# Every sacct query in bCheck is bounded by this timestamp. Job names
# embed $$ and pids recycle, so without it a name from an earlier run
# could be counted against this one.
#
export ATAC_START=$(date -d '-2 min' +%Y-%m-%dT%H:%M:%S)
RUNRE="^${TAG}_.*_$$\$"

#
# Run records, read by bin/checkRun.sh together with sacct.
#
#   00.RUNSTATUS.txt     KEY=VALUE. STATUS is RUNNING from here on. It is
#                        set to COMPLETED only as the last action of this
#                        script, after the deliverables check, and to
#                        FAILED (with the STAGE it failed in) by the EXIT
#                        trap. If the control job is killed outright
#                        (SIGKILL on its memory cap, node failure) the
#                        file stays at RUNNING and checkRun.sh finds the
#                        control job dead in sacct.
#
#   SLURM.CTRL/jobs.tsv  one line per stage job: job id, stage, sample
#                        and log path, written at submission.
#
RUNSTATUS=00.RUNSTATUS.txt
JOBS=SLURM.CTRL/jobs.tsv

CTRL_JOBID=${SLURM_JOB_ID:-none}
CTRL_LOG=none
if [ "$CTRL_JOBID" != "none" ]; then
    CTRL_LOG=$(scontrol show job $CTRL_JOBID 2>/dev/null | sed -n 's/^ *StdOut=//p')
    CTRL_LOG=${CTRL_LOG:-SLURM.CTRL/$CTRL_JOBID.out}
fi
echo CTRL_JOBID=$CTRL_JOBID

RUN_STATUS=RUNNING
RUN_STAGE=SETUP
RUN_RC=""
RUN_MESSAGE=""
RUN_STARTED=$(date '+%Y-%m-%d %H:%M:%S')
RUN_FINISHED=""
RUN_DONE=0
GENOME=""
SAMPLES=""

writeRunStatus() {
    cat >$RUNSTATUS.tmp <<EOF
STATUS=$RUN_STATUS
STAGE=$RUN_STAGE
RC=$RUN_RC
MESSAGE=$RUN_MESSAGE
STARTED=$RUN_STARTED
FINISHED=$RUN_FINISHED
VERSION=$SCRIPT_VERSION
CTRL_JOBID=$CTRL_JOBID
CTRL_LOG=$CTRL_LOG
HOST=$(hostname)
PID=$$
ATAC_START=$ATAC_START
GENOME=$GENOME
MAPQ=$MAPQ
SAMPLES=$SAMPLES
BAMS=$BAMS
JOBS=$JOBS
EOF
    mv $RUNSTATUS.tmp $RUNSTATUS
}

setStage() {
    RUN_STAGE=$1
    writeRunStatus
    echo
    echo "==== STAGE $RUN_STAGE $(date '+%Y-%m-%d %H:%M:%S')"
}

mkdir -p SLURM.CTRL
printf '#JOBID\tSTAGE\tSAMPLE\tLOG\n' >$JOBS
writeRunStatus

#
# Any exit before RUN_DONE=1 is a failure: an error under set -e, a failed
# bCheck, or scancel / time limit (SIGTERM). Record it in the status file
# and cancel the rest of this run's stage jobs so they do not keep running
# after the control job is gone. SIGTERM reaches the trap promptly only
# while the control job is in bSync (see slurmTools.sh), which is where it
# spends nearly all of its time; during the R reports at the end it may
# be SIGKILLed first, and then checkRun.sh finds the dead control job in
# sacct.
#
onSignal() {
    RUN_MESSAGE="control job got SIG$1 (scancel, time limit or node shutdown)"
    exit $2
}
trap 'onSignal TERM 143' TERM
trap 'onSignal INT 130' INT

onExit() {
    local rc=$?
    if [ -n "$BSYNC_PID" ]; then
        kill $BSYNC_PID 2>/dev/null || true
    fi
    if [ "$RUN_DONE" != "1" ]; then
        RUN_STATUS=FAILED
        RUN_RC=$rc
        RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
        writeRunStatus
        echo
        echo "    pipe.sh FAILED in stage $RUN_STAGE (rc=$rc); cancelling queued jobs of this run"
        echo
        bKill "$RUNRE"
    fi
}
trap onExit EXIT

GENOME=$($SDIR/bin/getGenomeBuildBAM.sh $1)

if [[ $GENOME =~ unknown ]]; then
    echo
    echo "    FATAL ERROR: UNKNOWN GENOME"
    echo "    "$GENOME
    echo
    exit 1
fi

echo GENOME=$GENOME

getSMTag () {
    samtools view -H $1 \
        | fgrep "@RG" \
        | head -1 \
        | tr '\t' '\n' \
        | fgrep SM: \
        | head -1 \
        | sed 's/SM://'
}

#
# The SM tag names the out/<SID>/ directory, the manifest lines and the
# per-sample deliverables, so a BAM without one, or two BAMs with the same
# one, cannot be processed.
#
for BAM in $BAMS; do
    SID=$(getSMTag $BAM)
    if [ -z "$SID" ]; then
        echo
        echo "    FATAL ERROR: no @RG SM tag in [$BAM]"
        echo
        exit 1
    fi
    SAMPLES="$SAMPLES $SID"
done
SAMPLES=${SAMPLES# }

DUPS=$(echo $SAMPLES | tr ' ' '\n' | sort | uniq -d)
if [ -n "$DUPS" ]; then
    echo
    echo "    FATAL ERROR: SM tag used by more than one BAM: "$DUPS
    echo
    exit 1
fi

echo SAMPLES=$SAMPLES
writeRunStatus

#
# Scrub the Slurm variables of the control job itself. bsub submits with
# the default --export=ALL, so they would otherwise ride along into every
# stage job. SLURM_CONF stays; the client tools need it.
#
for V in ${!SLURM_@}; do
    if [ "$V" != "SLURM_CONF" ]; then
        unset $V
    fi
done

#
# Short partitions for every stage: ~/bin/bsub sends anything under two
# hours to cmobic_short,cpushort. -M is the total memory for the job and
# a hard cgroup limit, unlike the LSF rusage[] request, which was a
# per-slot scheduling hint.
#
RUNTIME="-W 1:59:00"

#
# atacSub STAGE SAMPLE LOGDIR "BSUB_OPTS" CMD [ARG ...]
#
# Submit one stage job as ${TAG}_STAGE_$$ through bin/runStage.sh, which
# appends #ATAC_EXIT=<rc> to the job log, and record it in $JOBS.
#
atacSub() {
    local stage=$1
    local sample=$2
    local logdir=$3
    local opts=$4
    local out jobid
    shift 4

    if ! out=$(bsub $RUNTIME -o $logdir/ -J ${TAG}_${stage}_$$ $opts \
                   $SDIR/bin/runStage.sh "$@"); then
        echo "$out"
        echo
        echo "    FATAL: bsub failed for stage $stage sample $sample"
        echo
        exit 1
    fi
    echo "$out"

    jobid=$(echo "$out" | sed -n 's/^Job <\([0-9][0-9]*\)> is submitted.*/\1/p')
    if [ -z "$jobid" ]; then
        echo
        echo "    FATAL: no job id in bsub output for stage $stage sample $sample"
        echo
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\n' "$jobid" "$stage" "$sample" "$logdir/$jobid.out" >>$JOBS
}

#
# waitStage STAGE
#
# Block until every ${TAG}_STAGE_$$ job is done, then abort the run if any
# of them did not reach COMPLETED.
#
waitStage() {
    setStage $1
    bSync ${TAG}_$1_$$
    bCheck ${TAG}_$1_$$
}

setStage SUBMIT

for BAM in $BAMS; do
    atacSub POST2 $(getSMTag $BAM) SLURM.01.POST "-M 32G" \
        $SDIR/postMapBamProcessing_ATACSeq.sh -q $MAPQ $GENOME $BAM
done

waitStage POST2

for BEDZ in out/*/*.bed.gz; do
    SID=$(basename $(dirname $BEDZ))
    atacSub BW2 $SID SLURM.02.BW "-M 24G" \
        $SDIR/makeBigWigFromBEDZ.sh $GENOME $BEDZ
    atacSub CALLP2 $SID SLURM.03.CALLP "-n 3 -M 18G" \
        $SDIR/callPeaks_ATACSeq.sh $GENOME $BEDZ
done

for PBAM in out/*/*_postProcess.bam; do
    atacSub Index $(basename $(dirname $PBAM)) SLURM.04c.INDEX "-M 4G" \
        samtools index $PBAM
done

waitStage BW2
waitStage CALLP2

#
# MergePeaks -> Count -> DESEQ run in sequence through bSync rather than
# with -w post_done(): the bsub shim maps that to afterok but Slurm needs
# numeric job ids, and an afterok whose parent failed pends forever in
# DependencyNeverSatisfied.
#
atacSub MergePeaks all SLURM.04a.CALLP "-n 3 -M 24G" \
    $SDIR/mergePeaksToSAF.sh callpeaks \>macsPeaksMerged.saf

waitStage MergePeaks

PBAMS=$(ls out/*/*_postProcess.bam)
atacSub Count all SLURM.04b.CALLP "-n 10 -M 24G" \
    $SDIR/bin/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a macsPeaksMerged.saf \
        -o peaks_raw_fcCounts.txt \
        $PBAMS

waitStage Count

atacSub DESEQ all SLURM.05.DESEQ "-M 24G" \
    Rscript --no-save $SDIR/R/getDESeqScaleFactors.R

setStage MANIFEST

if [ ! -e "sampleManifest.csv" ]; then
    echo "MapID,SampleID,Group" > sampleManifest.csv
    for file in out/*/*bam; do getSMTag $file; done | sort >mapid
    cat mapid | sed 's/^s_//' >sid
    cat sid | perl -pe 's/(-|_)\d+$//' >gid
    paste mapid sid gid | tr '\t' ',' >> sampleManifest.csv
fi

waitStage DESEQ
waitStage Index

for PBAM in out/*/*_postProcess.bam; do
    atacSub TSSE $(basename $(dirname $PBAM)) SLURM.06.QC "-n 2 -M 32G" \
        $SDIR/bin/computeTSSEnrich.sh $PBAM
done

waitStage TSSE

setStage REPORTS

Rscript $SDIR/plotINSStats.R
Rscript $SDIR/R/analyzeATAC.R sampleManifest.csv

tee -a 00.POST_RUN.txt << 'EOF'

May want to check sampleManifest.csv
and rerun

    Rscript $SDIR/plotINSStats.R

    Rscript $SDIR/R/analyzeATAC.R sampleManifest.csv

and copy output to `atacSeq/metrics`

EOF

setStage STAGING

mkdir -p atacSeq/atlas
mkdir atacSeq/bigwig atacSeq/macs
mkdir -p atacSeq/metrics
mkdir -p out/postBams
mkdir out/metrics
mkdir out/bed

FAILED_JOBS=$(bCheckAll "$RUNRE")

if [ "$FAILED_JOBS" != "" ]; then
  echo -e "\n\n\nFailed Slurm jobs\n\n"
  echo "$FAILED_JOBS"
  echo -e "\n\n"
  exit 1
fi

mv macsPeaksMerged* atacSeq/atlas
cp peaks_raw_fcCounts.txt* atacSeq/atlas/
mv *_postProcess.shifted.10mNorm.bw atacSeq/bigwig
cp -val callpeaks/* atacSeq/macs
cp *__postInsDistribution.pdf *__ATACSeqQC.pdf atacSeq/metrics

cp -val out/*/*___INS.* atacSeq/metrics
cp -val out/*/*enrich* atacSeq/metrics

#
# Every sample must have each of its deliverables, not just the run as a
# whole. A peak file may legitimately be empty, so it only has to exist;
# everything else has to be non-empty.
#
setStage DELIVERABLES

MISSING=""

checkFile() {
    if [ ! -s "$1" ]; then
        MISSING="$MISSING $1"
    fi
}

checkFile atacSeq/atlas/macsPeaksMerged.saf
checkFile atacSeq/atlas/peaks_raw_fcCounts.txt
checkFile scaleFactorsDESeq2.csv

for PATTERN in "atacSeq/metrics/*__postInsDistribution.pdf" "atacSeq/metrics/*__ATACSeqQC.pdf"; do
    if ! compgen -G "$PATTERN" >/dev/null; then
        MISSING="$MISSING $PATTERN"
    fi
done

for SID in $SAMPLES; do
    checkFile atacSeq/bigwig/${SID}_postProcess.shifted.10mNorm.bw
    checkFile atacSeq/metrics/${SID}_postProcess___INS.txt
    checkFile atacSeq/metrics/${SID}_postProcess.tss_enrich.csv
    PEAKS=atacSeq/macs/${SID}_postProcess.shifted/${SID}_postProcess.shifted_peaks.narrowPeak
    if [ ! -e "$PEAKS" ]; then
        MISSING="$MISSING $PEAKS"
    fi
done

if [ -n "$MISSING" ]; then
    echo
    echo "    ERROR missing deliverables:"
    for F in $MISSING; do
        echo "        $F"
    done
    echo
    RUN_MESSAGE="missing deliverables:$MISSING"
    exit 1
fi

RUN_STATUS=COMPLETED
RUN_STAGE=DONE
RUN_RC=0
RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
writeRunStatus
RUN_DONE=1

echo
echo "==== ATAC-Seq run COMPLETED $RUN_FINISHED"
echo
