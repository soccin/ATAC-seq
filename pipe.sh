#!/bin/bash

# CMD:
#    bsub -n 1 -q control -o LSF.CTRL/ -J CTRL.ATAC ./pipe.sh
#

SDIR="$( cd "$( dirname "$0" )" && pwd )"

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

RUNTIME="-We 59"
echo $BAMS \
    | xargs -n 1 bsub $RUNTIME -o LSF.POST/ -J ${TAG}_POST2_$$ -R "rusage[mem=24]" -M 25 $SDIR/postMapBamProcessing_ATACSeq.sh

bSync ${TAG}_POST2_$$

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.BW/ -J ${TAG}_BW2_$$ $SDIR/makeBigWigFromBEDZ.sh
ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.CALLP/ -J ${TAG}_CALLP2_$$ -n 3 -R "rusage[mem=16]" -M 17 \
        $SDIR/callPeaks_ATACSeq.sh
