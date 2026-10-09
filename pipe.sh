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
    exit
}

if [ "$#" -lt "1" ]; then
    usage
fi

echo "MAPQ = ${MAPQ}"

BAMS=$*
echo SDIR=$SDIR
echo BAMS=$BAMS
echo VERSION=$SCRIPT_VERSION
echo CTRL_JOBID=${SLURM_JOB_ID:-none}

GENOME=$($SDIR/bin/getGenomeBuildBAM.sh $1)

if [[ $GENOME =~ unknown ]]; then
    echo
    echo "    FATAL ERROR: UNKNOWN GENOME"
    echo "    "$GENOME
    echo
    exit 1
fi

echo GENOME=$GENOME

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
# Every sacct query in bCheck is bounded by this timestamp. Job names
# embed $$ and pids recycle, so without it a name from an earlier run
# could be counted against this one.
#
export ATAC_START=$(date -d '-2 min' +%Y-%m-%dT%H:%M:%S)
RUNRE="^${TAG}_.*_$$\$"

#
# An aborted run must not leave the rest of a stage in the queue.
#
onExit() {
    local rc=$?
    if [ "$rc" != "0" ]; then
        echo
        echo "    pipe.sh FAILED (rc=$rc); cancelling queued jobs of this run"
        echo
        bKill "$RUNRE"
    fi
}
trap onExit EXIT

#
# Short partitions for every stage: ~/bin/bsub sends anything under two
# hours to cmobic_short,cpushort. -M is the total memory for the job and
# a hard cgroup limit, unlike the LSF rusage[] request, which was a
# per-slot scheduling hint.
#
RUNTIME="-W 1:59:00"

for BAM in $BAMS; do
    bsub $RUNTIME -o SLURM.01.POST/ -J ${TAG}_POST2_$$ -M 32G \
        $SDIR/postMapBamProcessing_ATACSeq.sh -q $MAPQ $GENOME $BAM
done

bSync ${TAG}_POST2_$$
bCheck ${TAG}_POST2_$$

for BEDZ in out/*/*.bed.gz; do
    bsub $RUNTIME -o SLURM.02.BW/ -J ${TAG}_BW2_$$ -M 24G \
        $SDIR/makeBigWigFromBEDZ.sh $GENOME $BEDZ
    bsub $RUNTIME -o SLURM.03.CALLP/ -J ${TAG}_CALLP2_$$ -n 3 -M 18G \
        $SDIR/callPeaks_ATACSeq.sh $GENOME $BEDZ
done

for PBAM in out/*/*_postProcess.bam; do
    bsub $RUNTIME -o SLURM.04c.INDEX/ -J ${TAG}_Index_$$ -M 4G \
        samtools index $PBAM
done

bSync ${TAG}_BW2_$$
bCheck ${TAG}_BW2_$$
bSync ${TAG}_CALLP2_$$
bCheck ${TAG}_CALLP2_$$

#
# MergePeaks -> Count -> DESEQ run in sequence through bSync rather than
# with -w post_done(): the bsub shim maps that to afterok but Slurm needs
# numeric job ids, and an afterok whose parent failed pends forever in
# DependencyNeverSatisfied.
#
bsub $RUNTIME -o SLURM.04a.CALLP/ -J ${TAG}_MergePeaks_$$ -n 3 -M 24G \
    $SDIR/mergePeaksToSAF.sh callpeaks \>macsPeaksMerged.saf

bSync ${TAG}_MergePeaks_$$
bCheck ${TAG}_MergePeaks_$$

PBAMS=$(ls out/*/*_postProcess.bam)
bsub $RUNTIME -o SLURM.04b.CALLP/ -J ${TAG}_Count_$$ -n 10 -M 24G \
    $SDIR/bin/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a macsPeaksMerged.saf \
        -o peaks_raw_fcCounts.txt \
        $PBAMS

bSync ${TAG}_Count_$$
bCheck ${TAG}_Count_$$

bsub $RUNTIME -o SLURM.05.DESEQ/ -J ${TAG}_DESEQ_$$ -M 24G \
    Rscript --no-save $SDIR/R/getDESeqScaleFactors.R

getSMTag () {
    samtools view -H $1 \
        | fgrep "@RG" \
        | head -1 \
        | tr '\t' '\n' \
        | fgrep SM: \
        | head -1 \
        | sed 's/SM://'
}


if [ ! -e "sampleManifest.csv" ]; then
    echo "MapID,SampleID,Group" > sampleManifest.csv
    for file in out/*/*bam; do getSMTag $file; done | sort >mapid
    cat mapid | sed 's/^s_//' >sid
    cat sid | perl -pe 's/(-|_)\d+$//' >gid
    paste mapid sid gid | tr '\t' ',' >> sampleManifest.csv
fi

bSync ${TAG}_DESEQ_$$
bCheck ${TAG}_DESEQ_$$

bSync ${TAG}_Index_$$
bCheck ${TAG}_Index_$$

for PBAM in out/*/*_postProcess.bam; do
    bsub $RUNTIME -o SLURM.06.QC/ -J ${TAG}_TSSE_$$ -n 2 -M 32G \
        $SDIR/bin/computeTSSEnrich.sh $PBAM
done

bSync ${TAG}_TSSE_$$
bCheck ${TAG}_TSSE_$$

Rscript $SDIR/plotINSStats.R
Rscript $SDIR/R/analyzeATAC.R sampleManifest.csv

tee -a 00.POST_RUN.txt << 'EOF'

May want to check sampleManifest.csv
and rerun

    Rscript $SDIR/plotINSStats.R

    Rscript $SDIR/R/analyzeATAC.R sampleManifest.csv

and copy output to `atacSeq/metrics`

EOF

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

if [ ! -e atacSeq/atlas/macsPeaksMerged.saf ]; then
    echo
    echo ERROR Postprocessing failed
    echo
    exit 1
fi
