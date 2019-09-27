#!/bin/bash

RESDIR=$1

echo $RESDIR

mkdir -p $RESDIR/atacSeq/atlas
mkdir -p $RESDIR/atacSeq/bigwig
mkdir -p $RESDIR/atacSeq/macs

rsync -rvP *.bw $RESDIR/atacSeq/bigwig
cp macsPeaksMerged.saf $RESDIR/atacSeq/atlas
cat macsPeaksMerged.saf | awk '{print $2,$3,$4,$1,".",$5}' | tr ' ' '\t' >$RESDIR/atacSeq/atlas/macsPeaksMerged.bed
rsync -rvP callpeaks/* $RESDIR/atacSeq/macs


#mkdir -p $RESDIR/atacSeq/diffpeaks
