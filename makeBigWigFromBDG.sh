#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

BDG=$1
OUT=$(basename $BDG | sed 's/.bdg/.bw/')

GENOME=$SDIR/mouse_mm10.genome

TDIR=/scratch/socci
mkdir -p $TDIR
TMP=$(mktemp -p $TDIR)
echo $TMP, $OUT
bedtools intersect -a $BDG -b ${GENOME}.bed -wa | sort -k1,1V -k2,2n -S 16G >$TMP
$SDIR/bedGraphToBigWig $TMP $GENOME $OUT


