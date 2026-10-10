#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

MERGEID=$1
SAMPLES=$2

SAMPLES=$(echo $SAMPLES | tr ',' ' ')
echo $SAMPLES

zcat $SAMPLES | sort -S 24g -T "$TMPDIR" -k1,1V -k2,2n | gzip -c - >merge_${MERGEID}__postProcess.shiftedPE.bed.gz
