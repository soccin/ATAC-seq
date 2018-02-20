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

#TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
TDIR=_scratch_ATACSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

samtools view -f 2 -F 524 $IBAM -u >$TDIR/step1.bam
picardV2 SortSam I=$TDIR/step1.bam O=$TDIR/step2.bam SO=queryname MAX_RECORDS_IN_RAM=5000000
samtools view -h $TDIR/step2.bam \
    | $SDIR/assign_multimappers.py -k 60 --paired-end \
    | samtools fixmate -r /dev/stdin $TDIR/step3.bam

#samtools view -f 2 -F 1804 -u $TDIR/step3.bam >$TDIR/step4.bam
samtools view -f 2 -F 3852 -u $TDIR/step3.bam >$TDIR/step4.bam

picardV2 SortSam I=$TDIR/step4.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000

#
# NDS add, remove non-standard chromosomes
#
# This code from:
#    https://www.biostars.org/p/187204/


bedtools bamtobed -i $TDIR/step4.bam -bedpe \
    | awk '$1 !~ /_/{print $0}' \
    | awk -F'\t' 'BEGIN {OFS = FS}{if($9=="+"){print $1,$2+4,$6+4,$7,$8,$9}else if($9=="-"){print $1,$2-5,$6-5,$7,$8,$9}}' \
    | sort -S20g -k1,1V -k2,2n \
    | gzip -nc >${OBAM/.bam/.shiftedPE.bed}.gz

#rm -rf $TDIR
