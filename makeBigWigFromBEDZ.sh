#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

set -o pipefail

GENOMEBUILD=$1
BEDZ=$2

$SDIR/bin/genomeInfo.sh "$GENOMEBUILD" TAG >/dev/null || exit 1
GENOME=$SDIR/lib/genomes/$GENOMEBUILD.genome
if [ ! -s "$GENOME" ]; then
    echo "FATAL ERROR: no genome file $GENOME"
    exit 1
fi

echo "md5sum ${BEDZ}"
echo $(md5sum ${BEDZ})

checkMD5=$(md5sum -c ${BEDZ}.md5 | awk '{print $2}'); 
if [ "$checkMD5" != "OK" ]; then
    echo
    echo md5check fail
    echo sleeping for 300 secs 
    sleep 300
    checkMD5=$(md5sum -c ${BEDZ}.md5 | awk '{print $2}'); 
    if [ "$checkMD5" != "OK" ]; then 
        echo 
        echo Second md5check fail
        echo FATAL ERROR Stopping
        echo
        exit 1
    fi
fi


if [ "$#" == "3" ]; then
    scaleFactor=$2
    echo "$BEDZ sizeFactorNorm scaleFactor "$scaleFactor
    OUT=$(basename $BEDZ | sed 's/.bed.gz/.sizeFactorNorm.bw/')
else
    count=$(zcat $BEDZ | cut -f4 | sort -S20g -T "$TMPDIR" | uniq | wc -l) \
        || { echo "FATAL ERROR: read count of $BEDZ failed"; exit 1; }
    if [ "$count" -eq 0 ]; then
        echo "FATAL ERROR: no reads in $BEDZ"
        exit 1
    fi
    scaleFactor=$(bc -l <<< "10000000/$count")
    echo "$BEDZ 10mNorm scaleFactor "$scaleFactor
    OUT=$(basename $BEDZ | sed 's/.bed.gz/.10mNorm.bw/')
fi

GENOME_BED=$(mktemp -p "$TMPDIR" genome.XXXXXX.bed) \
    || { echo "FATAL ERROR: mktemp in $TMPDIR failed"; exit 1; }
awk 'BEGIN {OFS = "\t"} {print $1, 0, $2}' $GENOME >$GENOME_BED \
    || { echo "FATAL ERROR: cannot write $GENOME_BED"; exit 1; }

#
# N.B. -scale argument is multiplicative
#
#    scaledCounts = rawCounts * $scaleFactor
#

echo "scaleFactor=${scaleFactor}"

zcat $BEDZ \
    | bedtools slop -i - -g $GENOME -s -l 0 -r 0 \
    | bedtools intersect -nonamecheck -a - -b $GENOME_BED \
    | bedtools genomecov -i - -g $GENOME -bg -scale $scaleFactor \
    | $SDIR/bin/wigToBigWig stdin $GENOME $OUT \
    || { echo "FATAL ERROR: bigWig build for $BEDZ failed"; rm -f $GENOME_BED; exit 1; }

rm -f $GENOME_BED
