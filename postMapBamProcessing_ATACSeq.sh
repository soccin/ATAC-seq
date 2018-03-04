#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

IBAM=$1

if [ "$#" == "1" ]; then
    OBAM=${IBAM/.bam/_postProcess.bam}
    OBAM=$(basename $OBAM)
    echo $OBAM
else
    OBAM=$2
fi

TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
#TDIR=_scratch_ATACSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

# f 3 ==> paired, proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
samtools view -q 10 -f 3 -F 1804 $IBAM -u >$TDIR/step1.bam
picardV2 SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000

#
# Do Tn5 shift and remove non-standard chromosomes
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | egrep -v "chrUn|_random" \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shifted.bed.gz}

rm -rf $TDIR
