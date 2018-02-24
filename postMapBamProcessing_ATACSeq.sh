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

#samtools view -f 2 -F 524 $IBAM -u >$TDIR/step1.bam

# f 66 ==> proper pair, first in pair
# F 524 ==> unmapped, mate unmapped, read fails QC
samtools view -f 66 -F 524 $IBAM -u >$TDIR/step1.bam

picardV2 SortSam I=$TDIR/step1.bam O=$TDIR/step2.bam SO=queryname MAX_RECORDS_IN_RAM=5000000
samtools view -h $TDIR/step2.bam \
    | $SDIR/assign_multimappers.py -k 60 --paired-end \
    | samtools fixmate -r /dev/stdin $TDIR/step3.bam

#samtools view -q 20 -f 2 -F 1804 -u $TDIR/step3.bam >$TDIR/step4.bam
samtools view -q 20 -F 1796 -u $TDIR/step3.bam >$TDIR/step4.bam
picardV2 SortSam I=$TDIR/step4.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000

#
# NDS add, remove non-standard chromosomes
#

bedtools bamtobed -i $OBAM \
    | awk -F'\t' 'BEGIN{OFS="\t"}{$5="1000";print $0}' \
    | awk '$1 !~ /_/{print $0}' \
    | gzip -nc >$TDIR/step5.bed.gz \

#
# Do Tn5 shift and get rid of non-standard chromosomes

zcat $TDIR/step5.bed.gz \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shifted.bed.gz}

rm -rf $TDIR
