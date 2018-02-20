#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

IBAM=$1

if [ "$#" == "1" ]; then
    OBAM=${IBAM/.bam/_postProcess.bam}
    echo $OBAM
else
    OBAM=$2
fi

TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

samtools view -f 2 -F 524 $IBAM -u >$TDIR/step1.bam
picardV2 SortSam I=$TDIR/step1.bam O=$TDIR/step2.bam SO=queryname
samtools view -h $TDIR/step2.bam \
    | $SDIR/assign_multimappers.py -k 60 --paired-end \
    | samtools fixmate -r /dev/stdin $TDIR/step3.bam
samtools view -f 2 -F 1804 -u $TDIR/step3.bam >$TDIR/step4.bam
picardV2 SortSam I=$TDIR/step4.bam O=$OBAM SO=coordinate

echo
echo "27"
echo

#
# NDS add, remove non-standard chromosomes
#

bedtools bamtobed -i $OBAM -bedpe -mate 1 \
    | awk '$1 !~ /_/{print $0}' \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($9 == "+") {$2 = $2 + 4} else if ($9 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shiftedPE.bed}

#rm -rf $TDIR
