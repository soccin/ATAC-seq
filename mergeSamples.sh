#!/bin/bash

MERGEID=$1
SAMPLES=$2

SAMPLES=$(echo $SAMPLES | tr ',' ' ')
echo $SAMPLES

zcat $SAMPLES | sort -S 24g -k1,1V -k2,2n | gzip -c - >merge_${MERGEID}__postProcess.shiftedPE.bed.gz

