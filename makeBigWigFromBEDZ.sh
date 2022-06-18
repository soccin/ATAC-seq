#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

GENOMEBUILD=$1
BEDZ=$2

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
    count=$(zcat $BEDZ | cut -f4 | sort -S20g | uniq | wc -l)
    scaleFactor=$(bc -l <<< "10000000/$count")
    echo "$BEDZ 10mNorm scaleFactor "$scaleFactor
    OUT=$(basename $BEDZ | sed 's/.bed.gz/.10mNorm.bw/')
fi

case $GENOMEBUILD in

    b37)
    GENOME=$SDIR/lib/genomes/human_b37.genome
    ;;

    mm10)
    GENOME=$SDIR/lib/genomes/mouse_mm10.genome
    ;;

    sCer+sMik_IFO1815)
    GENOME=$SDIR/lib/genomes/sCer+sMik_IFO1815.genome
    ;;

    *)
    echo
    echo "    Unknown GENOMEBUILD [$GENOMEBUILD]"
    echo
    exit 1
    ;;

esac

# TDIR=/scratch/socci
# mkdir -p $TDIR
# TMP=$(mktemp -p $TDIR)

#
# N.B. -scale argument is multiplicative
#
#    scaledCounts = rawCounts * $scaleFactor
#

echo "scaleFactor=${scaleFactor}"

zcat $BEDZ \
    | bedtools slop -i - -g $GENOME -s -l 0 -r 0 \
    | egrep -v "chrUn|_random|_unplaced|GL|NC_|hs37d5" \
    | bedtools genomecov -i - -g $GENOME -bg -scale $scaleFactor \
    | $SDIR/bin/wigToBigWig stdin $GENOME $OUT
