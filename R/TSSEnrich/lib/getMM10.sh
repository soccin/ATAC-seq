#!/bin/bash

set -euo pipefail

# Write mm10 TSS enrichment refs into this directory:
#   mm10_tss.bed       ENCODE mm10 unique GENCODE vM7 TSS
#   mm10.chrom.sizes   contig sizes from the given mm10 BAM or .dict
#
# ENCODE's file (encode-pipeline-genome-data/mm10/ataqc/
# mm10_gencode_tss_unique.bed.gz) returned 404 on 2026-10-11, as did the
# whole bucket. It is rebuilt here from the GENCODE vM7 GTF with the
# kundajelab/ataqc recipe ENCODE used (README.md, "build a TSS BED from
# scratch") followed by sort -u, which gives ENCODE's file byte for
# byte; the md5 is checked. The recipe's sort comes first so that, of
# the 9 genes that share a TSS with another, sort -u keeps the same one
# as ENCODE did (the first by the whole line, not the first in the GTF).
#
# The header must be the one getGenomeBuildBAM.sh calls mm10 (IRIS
# M.musculus/mm10/mm10.dict, 66 contigs). chrom.sizes has to come from
# that header, not from ENCODE/UCSC, or TSS names will not join.

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/../../../bin/loadTools.sh

if [ "$#" != "1" ]; then
    echo "usage: getMM10.sh BAM|DICT" >&2
    exit 1
fi

HEADER=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
cd "$SDIR"

GTF_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M7/gencode.vM7.annotation.gtf.gz
GTF_GZ=$TMPDIR/gencode.vM7.annotation.gtf.gz
TSS_MD5=b0a24d5c8165f8782c5c6c6c98ab511f

curl -fsSL -o "$GTF_GZ" "$GTF_URL"

zcat "$GTF_GZ" \
    | grep -P '\tgene\t' \
    | grep 'protein_coding' \
    | grep -v 'level 3' \
    | awk -F '[\t|"]' '{ print $1"\t"$4"\t"$5"\t"$10"\t0\t"$7 }' \
    | awk -F '\t' 'BEGIN{ OFS="\t" } { if ($6=="+") { $3=$2-1; $2=$2-2 } else { $2=$3; $3=$3+1 } print }' \
    | sort -T "$TMPDIR" -k1,1 -k2,2n \
    | sort -T "$TMPDIR" -k1,1 -k2,2n -u \
    > mm10_tss.bed
rm -f "$GTF_GZ"

if [ "$(md5sum mm10_tss.bed | awk '{print $1}')" != "$TSS_MD5" ]; then
    echo "ERROR: mm10_tss.bed md5 is not $TSS_MD5" >&2
    exit 1
fi

case $HEADER in
    *.dict)
    cat "$HEADER"
    ;;
    *)
    samtools view -H "$HEADER"
    ;;
esac \
    | awk -F'\t' '$1=="@SQ" {
        sn=""; ln="";
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) sn=substr($i,4);
            if ($i ~ /^LN:/) ln=substr($i,4);
        }
        if (sn != "" && ln != "") print sn"\t"ln
    }' \
    | sort -T "$TMPDIR" -k2,2nr \
    > mm10.chrom.sizes
