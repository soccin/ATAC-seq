#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

set -o pipefail

if [ "$#" != "2" ]; then
    echo "usage: mergePeaksToSAF.sh GENOME MACS_OUTDIR >atacSeq/atlas/macsPeaksMerged.saf" >&2
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

#
# bedtools merge needs the peaks of all samples sorted as one stream.
# One sort reads every file, its names given on stdin (--files0-from),
# not as arguments: xargs would split a long file list into several
# sorts, each sorted on its own, and merge would stop at the first
# out-of-order peak.
#
# The merged peaks are BED (0-based start, end exclusive); SAF, as
# featureCounts reads it, is 1-based with both ends inclusive, so the
# start gets 1 added and the end stays. Until v2.0.0 the BED start was
# written as is, and featureCounts counted one base more on the left.
#
find $MACS_OUTDIR -name "$PATTERN" -print0 \
    | sort -S 16g -T "$TMPDIR" -k1,1V -k2,2n --files0-from=- \
    | bedtools merge -i - -d 500 \
    | bedtools subtract -nonamecheck -a - -b $BLACKLIST \
    | awk 'BEGIN {OFS = "\t"} {print "Peak_"++s, $1, $2 + 1, $3, "+"}'
