#!/bin/bash

# CMD:
#    bsub -n 1 -q control -o LSF.CTRL/ -J CTRL.ATAC ./pipe.sh
#

SDIR="$( cd "$( dirname "$0" )" && pwd )"

module load bedtools/2.27.1

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

GENOME=$($SDIR/getGenomeBuildBAM.sh $1)

if [[ $GENOME =~ unknown ]]; then
    echo
    echo "    FATAL ERROR: UNKNOWN GENOME"
    echo "    "$GENOME
    echo
    exit 1
fi

RUNTIME="-W 359"
RUNTIME_SHORT="-W 59"

echo $BAMS \
    | xargs -n 1 bsub $RUNTIME -o LSF.01.POST/ -J ${TAG}_POST2_$$ -R "rusage[mem=24]" \
        $SDIR/postMapBamProcessing_ATACSeq.sh -q $MAPQ

bSync ${TAG}_POST2_$$

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.02.BW/ -J ${TAG}_BW2_$$ -R "rusage[mem=24]" $SDIR/makeBigWigFromBEDZ.sh $GENOME

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME_SHORT -o LSF.03.CALLP/ -J ${TAG}_CALLP2_$$ -n 3 -R "rusage[mem=6]" \
        $SDIR/callPeaks_ATACSeq.sh $GENOME

bSync ${TAG}_CALLP2_$$

bsub $RUNTIME_SHORT -o LSF.04a.CALLP/ -J ${TAG}_MergePeaks_$$ -n 3 -R "rusage[mem=24]" \
    $SDIR/mergePeaksToSAF.sh callpeaks \>macsPeaksMerged.saf

PBAMS=$(ls *_postProcess.bam)
bsub $RUNTIME -o LSF.04b.CALLP/ -J ${TAG}_Count_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_MergePeaks_$$)" \
    $SDIR/bin/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a macsPeaksMerged.saf \
        -o peaks_raw_fcCounts.txt \
        $PBAMS

bsub $RUNTIME_SHORT -o LSF.05.DESEQ/ -J ${TAG}_DESEQ_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_Count_$$)" \
    Rscript --no-save $SDIR/R/getDESeqScaleFactors.R

echo "SampleID,Group,MapID" > sampleManifest.csv
ls *postProcess.bam | sed 's/_postProcess.bam//' | sed 's/.*_s_/s_/' | sort >mapid
cat mapid | sed 's/^s_//' >sid
cat sid | perl -pe 's/(-|_)\d+$//' >gid
paste sid gid mapid | tr '\t' ',' >> sampleManifest.csv

bSync ${TAG}_MergePeaks_$$
bSync ${TAG}_Count_$$
bSync ${TAG}_DESEQ_$$

Rscript --no-save $SDIR/plotINSStats.R
Rscript --no-save $SDIR/R/analyzeATAC.R sampleManifest.csv

mkdir -p atacSeq/atlas
mkdir atacSeq/bigwig atacSeq/macs
mkdir -p atacSeq/metrics
mkdir -p out/postBams
mkdir out/metrics
mkdir out/bed

# mv macsPeaksMerged* atacSeq/atlas
# mv *_postProcess.shifted.10mNorm.bw atacSeq/bigwig
# cp -val callpeaks/* atacSeq/macs
# mv *_postProcess.bam out/postBams
# mv *___INS.* out/metrics/
# mv *shifted.bed.gz out/bed
# mv *shifted.bed.gz.md5 out/bed
# mv *__postInsDistribution.pdf *__ATACSeqQC.pdf atacSeq/metrics

module unload bedtools
