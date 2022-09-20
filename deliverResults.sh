#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

RESDIR=$1

if [ "$#" != "1" ]; then
    echo
    echo usage: deliveryResults.sh RESDIR
    echo "    RESDIR=/ifs/res/seq/pi/invest/r_###"
    echo
    exit
fi

echo \$RESDIR=$RESDIR $(realpath $RESDIR)

echo
echo "sudo needed to set permissions"
echo
sudo chmod g+ws $RESDIR

if [ -e "macsPeaksMerged.saf" ]; then
    echo
    echo "Need to run postProcessing file move/copy"
    echo

    mv macsPeaksMerged* atacSeq/atlas
    mv *_postProcess.shifted.10mNorm.bw atacSeq/bigwig
    cp -val callpeaks/* atacSeq/macs
    cp *__postInsDistribution.pdf *__ATACSeqQC.pdf atacSeq/metrics

    cp -val out/*/*___INS.* atacSeq/metrics

fi

if [ ! -e atacSeq/atlas/macsPeaksMerged.saf ]; then
    echo
    echo ERROR Postprocessing failed
    echo
    exit 1
fi


rsync -avP atacSeq $RESDIR

ATAC_PROJECT_NUM=$($SDIR/extractProjectIDFromPath.py $(realpath $RESDIR) | sed 's/^Proj_//')

cat $SDIR/tmpldeliveryEmail.txt | sed 's@ATAC_PROJECT_NUM@'"$ATAC_PROJECT_NUM"'@'
