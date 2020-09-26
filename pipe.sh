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

TAG=q$PIPENAME

COMMAND_LINE=$*

function usage {
    echo
    echo "usage: $PIPENAME/pipe.sh BAM1 [BAM2 ... BAMN]"
    echo "version=$SCRIPT_VERSION"
    echo ""
    echo
    exit
}

if [ "$#" -lt "1" ]; then
    usage
fi

BAMS=$*
echo SDIR=$SDIR
echo BAMS=$BAMS

GENOME=$($SDIR/getGenomeBuildBAM.sh $1)

RUNTIME="-W 59"
echo $BAMS \
    | xargs -n 1 bsub $RUNTIME -o LSF.01.POST/ -J ${TAG}_POST2_$$ -R "rusage[mem=24]" $SDIR/postMapBamProcessing_ATACSeq.sh

bSync ${TAG}_POST2_$$

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.02.BW/ -J ${TAG}_BW2_$$ -R "rusage[mem=24]" $SDIR/makeBigWigFromBEDZ.sh $GENOME

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.03.CALLP/ -J ${TAG}_CALLP2_$$ -n 3 -R "rusage[mem=24]" \
        $SDIR/callPeaks_ATACSeq.sh

bSync ${TAG}_CALLP2_$$

bsub $RUNTIME -o LSF.04a.CALLP/ -J ${TAG}_MergePeaks_$$ -n 3 -R "rusage[mem=24]" \
    $SDIR/mergePeaksToSAF.sh callpeaks \>macsPeaksMerged.saf

PBAMS=$(ls *_postProcess.bam)
bsub $RUNTIME -o LSF.04b.CALLP/ -J ${TAG}_Count_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_MergePeaks_$$)" \
    $SDIR/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a macsPeaksMerged.saf \
        -o peaks_raw_fcCounts.txt \
        $PBAMS

bsub $RUNTIME -o LSF.05.DESEQ/ -J ${TAG}_DESEQ_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_Count_$$)" \
    Rscript --no-save $SDIR/getDESeqScaleFactors.R

module unload bedtools
