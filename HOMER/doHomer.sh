#!/bin/bash

DPEAKS=$1
COMPARISON=$2

HOMERFILE=$(Rscript --no-save getHomerBED.R $DPEAKS $((COMPARISON+1)))

echo $HOMERFILE

PATH=$PATH:/home/socci/Work/Users/SchmittA/MitobeY/B-101-131/Homer/bin/

mkdir -p outHomer

findMotifsGenome.pl $HOMERFILE hg19 outHomer/$(basename ${HOMERFILE/.tsv/})

