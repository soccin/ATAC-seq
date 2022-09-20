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

# mv macsPeaksMerged* atacSeq/atlas
# mv *_postProcess.shifted.10mNorm.bw atacSeq/bigwig
# cp -val callpeaks/* atacSeq/macs
# cp *__postInsDistribution.pdf *__ATACSeqQC.pdf atacSeq/metrics

# mv *_postProcess.bam out/postBams
# mv *___INS.* out/metrics/
# mv *shifted.bed.gz out/bed
# mv *shifted.bed.gz.md5 out/bed

rsync -avP atacSeq $RESDIR

ATAC_PROJECT_NUM=$($SDIR/extractProjectIDFromPath.py $(realpath $RESDIR) | sed 's/^Proj_//')

cat $SDIR/tmpldeliveryEmail.txt | sed 's@ATAC_PROJECT_NUM@'"$ATAC_PROJECT_NUM"'@'
