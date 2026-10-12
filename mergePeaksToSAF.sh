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
# The peaks are already blacklist-filtered (callPeaks_ATACSeq.sh), but
# merging within 500 bp can bridge a small blacklisted region that lies
# between two kept peaks; subtract the blacklist again so no atlas peak
# covers any of it.
#
find $MACS_OUTDIR \
    | fgrep narrowPeak \
    | xargs sort -S 16g -T "$TMPDIR" -k1,1V -k2,2n \
    | bedtools merge -i - -d 500 \
    | bedtools subtract -nonamecheck -a - -b $BLACKLIST \
    | awk '{print "Peak_"++s,$0,"+"}' | tr ' ' '\t'
