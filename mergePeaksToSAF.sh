#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

set -o pipefail

if [ "$#" != "2" ]; then
    echo "usage: mergePeaksToSAF.sh GENOME MACS_OUTDIR >macsPeaksMerged.saf" >&2
    exit 1
fi

GENOMEBUILD=$1
MACS_OUTDIR=$2

BLACKLIST=$SDIR/lib/genomes/$GENOMEBUILD.blacklist.bed
if [ ! -s "$BLACKLIST" ]; then
    echo "FATAL ERROR: no blacklist $BLACKLIST" >&2
    exit 1
fi

#
# Only the blacklist-filtered peaks of callPeaks_ATACSeq.sh
# (*_peaks.blacklistFiltered.narrowPeak), never MACS2's own
# *_peaks.narrowPeak in the same directory. Merging within 500 bp can
# still bridge a small blacklisted region that lies between two kept
# peaks; subtract the blacklist again so no atlas peak covers any of it.
#
PATTERN='*_peaks.blacklistFiltered.narrowPeak'

if [ -z "$(find $MACS_OUTDIR -name "$PATTERN")" ]; then
    echo "FATAL ERROR: no $PATTERN under $MACS_OUTDIR" >&2
    exit 1
fi

find $MACS_OUTDIR -name "$PATTERN" \
    | xargs sort -S 16g -T "$TMPDIR" -k1,1V -k2,2n \
    | bedtools merge -i - -d 500 \
    | bedtools subtract -nonamecheck -a - -b $BLACKLIST \
    | awk '{print "Peak_"++s,$0,"+"}' | tr ' ' '\t'
