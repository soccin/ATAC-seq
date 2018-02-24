#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

BEDZ=$1
OUT=$(basename $BEDZ | sed 's/.bed.gz/.bw/')

GENOME=$SDIR/mouse_mm10.genome

TDIR=/scratch/socci
mkdir -p $TDIR
TMP=$(mktemp -p $TDIR)
echo $TMP, $OUT
count=$(zcat $BEDZ | wc -l)

zcat $BEDZ \
    | bedtools genomecov -bg -g $GENOME -i - \
    | awk '{print $1,$2,$3,10000000*$4/'$count'}' | tr ' ' '\t' \
    >${TMP}

$SDIR/bedGraphToBigWig $TMP $GENOME $OUT

