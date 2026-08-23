#!/bin/bash

set -euo pipefail

# Write b37 TSS enrichment refs into this directory:
#   b37_tss.bed       ENCODE hg19 unique GENCODE TSS, chr prefix stripped
#   b37.chrom.sizes   contig sizes from the given BAM header
#
# The BAM must be a b37 alignment (contigs 1, 2, ..., MT). chrom.sizes has to
# come from that header, not from ENCODE/UCSC, or TSS names will not join.

SDIR=$(cd "$(dirname "$0")" && pwd)
cd "$SDIR"

if [ "$#" != "1" ]; then
    echo "usage: getB37.sh BAM" >&2
    exit 1
fi

BAM=$1
TSS_URL=https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz
TSS_GZ=hg19_gencode_tss_unique.bed.gz

curl -fsSL -o "$TSS_GZ" "$TSS_URL"

zcat "$TSS_GZ" \
    | sed 's/^chr//' \
    | sort -k1,1V -k2,2n \
    > b37_tss.bed
rm -f "$TSS_GZ"

samtools view -H "$BAM" \
    | awk -F'\t' '$1=="@SQ" {
        sn=""; ln="";
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) sn=substr($i,4);
            if ($i ~ /^LN:/) ln=substr($i,4);
        }
        if (sn != "" && ln != "") print sn"\t"ln
    }' \
    > b37.chrom.sizes
