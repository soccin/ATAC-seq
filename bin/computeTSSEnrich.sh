#!/bin/bash

set -euo pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)
RDIR=$SDIR/../R

BAM=$1
ODIR=$(dirname $BAM)

GENOME=$($SDIR/getGenomeBuildBAM.sh $BAM)
echo "\$GENOME=$GENOME"

READ_LEN=$(
    samtools view $BAM \
        | cut -f6 \
        | grep -P "^[0-9]+M$" \
        | tr -d "M" \
        | head -100000 \
        | sort -V \
        | uniq -c \
        | sort -nr \
        | head -1 \
        | awk '{print $2}'
)
echo "\$READ_LEN=$READ_LEN"

Rscript $RDIR/TSSEnrich/tss_enrich.R \
    --nodup-bam $BAM \
    --read-len $READ_LEN \
    --chrsz $RDIR/TSSEnrich/lib/${GENOME}.chrom.sizes \
    --tss $RDIR/TSSEnrich/lib/${GENOME}_tss.bed \
    --out-dir $ODIR
