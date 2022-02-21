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

echo $RESDIR

echo
echo "sudo needed to set permissions"
echo
sudo chmod g+ws $RESDIR

mkdir -p $RESDIR/atacSeq/atlas
mkdir -p $RESDIR/atacSeq/bigwig
mkdir -p $RESDIR/atacSeq/macs

rsync -rvP *.bw $RESDIR/atacSeq/bigwig
cp macsPeaksMerged.saf $RESDIR/atacSeq/atlas
cat macsPeaksMerged.saf | awk '{print $2,$3,$4,$1,".",$5}' | tr ' ' '\t' >$RESDIR/atacSeq/atlas/macsPeaksMerged.bed
rsync -rvP callpeaks/* $RESDIR/atacSeq/macs


#mkdir -p $RESDIR/atacSeq/diffpeaks

ATAC_PROJECT_NUM=$($SDIR/extractProjectIDFromPath.py $RESDIR | sed 's/^Proj_//')

cat $SDIR/tmpldeliveryEmail.txt | sed 's@ATAC_PROJECT_NUM@'"$ATAC_PROJECT_NUM"'@'
