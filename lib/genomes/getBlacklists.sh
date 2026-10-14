#!/bin/bash

set -euo pipefail

# Write the peak blacklists into this directory, the files the ENCODE
# ATAC-seq pipeline uses (scripts/download_genome_data.sh, VER=v3):
#   hg38.blacklist.bed   ENCFF356LFX (Boyle lab blacklist v2), as is
#   mm10.blacklist.bed   ENCFF547MET (Boyle lab blacklist v2), as is
#   b37.blacklist.bed    ENCFF001TDO, the hg19 list ENCODE uses
#                        (wgEncodeDacMapabilityConsensusExcludable), with
#                        the chr prefix stripped and chrM dropped
#
# GRCh37 and hg19 have the same coordinates on 1-22, X and Y. Their
# mitochondrial sequences differ (hg19 chrM is not the b37 MT), and
# mitochondrial reads are not in peaks anyway, so chrM is dropped.
#
# The md5 of each download is the one on the ENCODE portal
# (https://www.encodeproject.org/files/<ACCESSION>/, checked 2026-10-12).
# Every contig of each blacklist must be in lib/genomes/<TAG>.genome.
#
# All three files are built and checked in a temp directory first. The
# files in this directory are replaced only if every download, md5 and
# contig check passes; on any failure they are left as they were.

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/../../bin/loadTools.sh

cd "$SDIR"

WORK=$(mktemp -d "$TMPDIR/getBlacklists.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

URL=https://www.encodeproject.org/files

get() {
    local ACC=$1
    local MD5=$2
    local GZ=$WORK/$ACC.bed.gz
    curl -fsSL -o "$GZ" "$URL/$ACC/@@download/$ACC.bed.gz"
    if [ "$(md5sum "$GZ" | awk '{print $1}')" != "$MD5" ]; then
        echo "ERROR: $ACC.bed.gz md5 is not $MD5" >&2
        exit 1
    fi
    zcat "$GZ"
}

get ENCFF356LFX 393688b4f06c9ce26165d47433dd8c37 >$WORK/hg38.blacklist.bed
get ENCFF547MET 5d19eeb4197ae10540252c61399f3131 >$WORK/mm10.blacklist.bed
get ENCFF001TDO 895f1dd97b001a338eb25d665c23d571 \
    | sed 's/^chr//' \
    | awk -F'\t' '$1 != "M"' \
    >$WORK/b37.blacklist.bed

for TAG in b37 hg38 mm10; do
    if [ ! -s $WORK/$TAG.blacklist.bed ]; then
        echo "ERROR: $TAG.blacklist.bed is empty" >&2
        exit 1
    fi
    BAD=$(cut -f1 $WORK/$TAG.blacklist.bed | sort -u -T "$TMPDIR" \
        | awk 'NR == FNR {ok[$1] = 1; next} !($1 in ok)' $TAG.genome - \
        | tr '\n' ' ')
    if [ -n "$BAD" ]; then
        echo "ERROR: $TAG.blacklist.bed contigs not in $TAG.genome: $BAD" >&2
        exit 1
    fi
done

#
# Every check passed. Copy each file in beside its target, then rename
# it, so a target is never partly written.
#

for TAG in b37 hg38 mm10; do
    cp $WORK/$TAG.blacklist.bed $TAG.blacklist.bed.new
    mv $TAG.blacklist.bed.new $TAG.blacklist.bed
    echo "$TAG.blacklist.bed $(wc -l <$TAG.blacklist.bed) regions"
done
