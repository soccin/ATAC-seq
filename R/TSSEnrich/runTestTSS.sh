#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo ""
    echo "   Usage: runTestTSS.sh COORD_SORTED_BAM"
    echo ""
    exit 1
fi

BAM=$1
NAMESORT_BAM="${BAM%.bam}_namesort.bam"
SCRIPT_DIR="$(dirname "$0")"

echo "Input BAM:      $BAM"
echo "Name-sorted:    $NAMESORT_BAM"

samtools sort -n -@ 4 "$BAM" -o "$NAMESORT_BAM"
echo "Name-sort done."

Rscript "$SCRIPT_DIR/testTSS.R" "$NAMESORT_BAM"

rm "$NAMESORT_BAM"
echo "Cleaned up $NAMESORT_BAM"
