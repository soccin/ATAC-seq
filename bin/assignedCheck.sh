#!/bin/bash
#
# assignedCheck.sh BAM SAF OUT
#
# Check B, one sample: count the fragments of BAM that overlap any peak
# of the atlas SAF by at least one base, without featureCounts, and write
#
#   BAM<TAB>FRAGMENTS
#
# to OUT. bin/checkCounts.sh compares FRAGMENTS with the featureCounts
# Assigned total of the same BAM; they must be equal.
#
# A fragment is counted once if either of its reads overlaps a peak, as
# featureCounts -p --countReadPairs -O counts it: the number of distinct
# read names among the reads samtools view -L reports. samtools tests the
# whole reference span of a read, featureCounts its aligned blocks; they
# differ only for a deletion that covers a whole peak, which a short read
# cannot have.
#
# BAM is a postProcess BAM: already filtered to MAPQ, proper pairs and no
# secondary, duplicate or QC-fail reads, so no filter is applied here. If
# that filter ever changes, apply the same one here.
#
# The SAF is 1-based and inclusive; the BED made from it starts one base
# lower, so both tools test the same interval.
#

set -o pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/loadTools.sh

if [ "$#" != "3" ]; then
    echo
    echo "    usage: assignedCheck.sh BAM SAF OUT"
    echo
    exit 1
fi

BAM=$1
SAF=$2
OUT=$3

BED=$(mktemp -p "$TMPDIR" atlas.XXXXXX.bed) || exit 1

awk 'BEGIN {OFS = "\t"} {s = $3 - 1; if (s < 0) s = 0; print $2, s, $4, $1}' $SAF >$BED \
    || { echo "FATAL ERROR: cannot convert $SAF to BED"; exit 1; }

NFRAG=$(
    samtools view -L $BED $BAM \
        | cut -f1 \
        | sort -u -S 4G --parallel=${SLURM_CPUS_PER_TASK:-1} -T "$TMPDIR" \
        | wc -l
) || { echo "FATAL ERROR: fragment count of $BAM failed"; exit 1; }

rm -f $BED

printf '%s\t%s\n' "$BAM" "$NFRAG" >$OUT.tmp \
    && mv $OUT.tmp $OUT \
    || { echo "FATAL ERROR: cannot write $OUT"; exit 1; }

echo "assignedCheck: $BAM fragments over the atlas: $NFRAG"
