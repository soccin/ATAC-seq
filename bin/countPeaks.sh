#!/bin/bash
#
# countPeaks.sh MAPQ SAF OUT BAM [BAM ...]
#
# Count stage: fragments per atlas peak with the vendored featureCounts
# (subread v2.1.1), then Check A (bin/checkCounts.sh) on the result.
# Exits nonzero if either fails, so the run stops at STAGE=Count.
#
#   -p --countReadPairs  count fragments. In 2.x -p alone only says the
#                        input is paired, and every read is counted:
#                        every column doubles, with no error.
#   -O                   a fragment over two peaks counts in both
#   -Q MAPQ              the MAPQ of postMapBamProcessing_ATACSeq.sh
#   --tmpDir "$TMPDIR"   node-local per-job directory (bin/runStage.sh)
#
# Defaults that matter: singletons count as fragments (no -B), chimeric
# pairs count (no -C), multimappers do not.
#
# featureCounts v1.6.0, used until v2.0.0, sometimes did not pair the
# reads of one input file when run with several threads on several BAMs,
# and counted that sample about twice, exiting 0. Check A catches that.
#
# Run from the analysis directory; BAM paths become the column names.
#

set -o pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/loadTools.sh

if [ "$#" -lt "4" ]; then
    echo
    echo "    usage: countPeaks.sh MAPQ SAF OUT BAM [BAM ...]"
    echo
    exit 1
fi

MAPQ=$1
SAF=$2
OUT=$3
shift 3

THREADS=${SLURM_CPUS_PER_TASK:-1}

$SDIR/featureCounts -v

$SDIR/featureCounts -O -Q $MAPQ -p --countReadPairs -T $THREADS \
    --tmpDir "$TMPDIR" \
    -F SAF -a $SAF \
    -o $OUT \
    "$@" \
    || { echo "FATAL ERROR: featureCounts failed"; exit 1; }

$SDIR/checkCounts.sh $OUT.summary \
    || { echo "FATAL ERROR: Check A failed on $OUT.summary"; exit 1; }
