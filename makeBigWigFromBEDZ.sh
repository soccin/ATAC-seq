#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

BEDZ=$1

if [ "$#" == "2" ]; then
    scaleFactor=$2
    echo "$BEDZ sizeFactorNorm scaleFactor "$scaleFactor
    OUT=$(basename $BEDZ | sed 's/.bed.gz/.sizeFactorNorm.bw/')
else
    count=$(zcat $BEDZ | cut -f4 | sort -S20g | uniq | wc -l)
    scaleFactor=$(bc -l <<< "10000000/$count")
    echo "$BEDZ 10mNorm scaleFactor "$scaleFactor
    OUT=$(basename $BEDZ | sed 's/.bed.gz/.10mNorm.bw/')
fi

GENOME=$SDIR/mouse_mm10.genome

# TDIR=/scratch/socci
# mkdir -p $TDIR
# TMP=$(mktemp -p $TDIR)

#
# N.B. -scale argument is multiplicative
#
#    scaledCounts = rawCounts * $scaleFactor
#

zcat $BEDZ \
    | bedtools slop -i - -g $GENOME -s -l 0 -r 0 \
    | egrep -v "chrUn|_random" \
    | bedtools genomecov -i - -g $GENOME -bg -scale $scaleFactor \
    | $SDIR/wigToBigWig stdin $GENOME $OUT

