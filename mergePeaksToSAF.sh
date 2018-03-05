#!/bin/bash

MACS_OUTDIR=$1

find $MACS_OUTDIR \
    | fgrep narrowPeak \
    | xargs sort -S 16g -k1,1V -k2,2n \
    | bedtools merge -i - -d 500 \
    | awk '{print "Peak_"++s,$0,"+"}' | tr ' ' '\t'

