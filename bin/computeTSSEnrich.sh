#!/bin/bash

set -euo pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)
RDIR=$SDIR/../R

. $SDIR/loadTools.sh

#
# GENOME is the run's genome, detected once by pipe.sh from the input
# BAMs, not detected again here.
#
if [ "$#" != "2" ]; then
    echo "usage: computeTSSEnrich.sh GENOME BAM" >&2
    exit 1
fi

GENOME=$1
BAM=$2
ODIR=$(dirname $BAM)

$SDIR/genomeInfo.sh "$GENOME" TAG >/dev/null
echo "\$GENOME=$GENOME"

CHRSZ=$RDIR/TSSEnrich/lib/${GENOME}.chrom.sizes
TSS=$RDIR/TSSEnrich/lib/${GENOME}_tss.bed
for F in $CHRSZ $TSS; do
    if [ ! -s "$F" ]; then
        echo "ERROR: no TSS enrichment file $F" >&2
        exit 1
    fi
done

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
        | awk '{print $2}' \
    || true
)
echo "\$READ_LEN=$READ_LEN"

if [[ -z $READ_LEN ]]; then
    echo "ERROR: could not determine read length from BAM: $BAM" >&2
    exit 1
fi

Rscript $RDIR/TSSEnrich/tss_enrich.R \
    --nodup-bam $BAM \
    --read-len $READ_LEN \
    --chrsz $CHRSZ \
    --tss $TSS \
    --out-dir $ODIR
