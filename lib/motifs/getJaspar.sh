#!/bin/bash

set -euo pipefail

# Write the motif collection of the motif enrichment step into this
# directory:
#   JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
#       JASPAR 2024 CORE, vertebrates, non-redundant: 879 position
#       frequency matrices in JASPAR format, as is
#
# One collection serves human and mouse: JASPAR CORE vertebrate
# profiles are not species-specific.
#
# JASPAR publishes no md5. The one below is of the file downloaded on
# 2026-10-23; a different md5 means JASPAR changed the file, and it must
# be looked at before the pinned value is changed.
#
# The file is downloaded and checked in a temp directory first. The
# file in this directory is replaced only if the download, md5 and
# profile count pass; on any failure it is left as it was.

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/../../bin/loadTools.sh

cd "$SDIR"

WORK=$(mktemp -d "$TMPDIR/getJaspar.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

FILE=JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
URL=https://jaspar.elixir.no/download/data/2024/CORE/$FILE
MD5=60d87dabb351150bd6394f35874ca2ea
NPROFILES=879

curl -fsSL -o $WORK/$FILE "$URL"

if [ "$(md5sum $WORK/$FILE | awk '{print $1}')" != "$MD5" ]; then
    echo "ERROR: $FILE md5 is not $MD5" >&2
    exit 1
fi

N=$(grep -c "^>" $WORK/$FILE)
if [ "$N" != "$NPROFILES" ]; then
    echo "ERROR: $FILE has $N profiles, not $NPROFILES" >&2
    exit 1
fi

#
# Every check passed. Copy the file in beside its target, then rename
# it, so the target is never partly written.
#

cp $WORK/$FILE $FILE.new
mv $FILE.new $FILE
echo "$FILE $N profiles"
