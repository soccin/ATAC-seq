#!/bin/bash
#
# loadTools.sh -- put the external tools the pipeline needs on PATH
#
# Usage: . $SDIR/bin/loadTools.sh
#
# samtools comes from Lmod. Batch jobs inherit the exported module
# function through --export=ALL, but a bare shell may not have it, so
# source the init script when it is missing.
#
# bedtools and picard are expected on PATH (~/bin on IRIS): there is no
# bedtools module and no system picard wrapper. featureCounts and
# wigToBigWig are vendored in bin/ and called by explicit path.
#

ATAC_SAMTOOLS_MODULE=${ATAC_SAMTOOLS_MODULE:-samtools/1.20}

if [ "$(type -t module)" != "function" ]; then
    . /etc/profile.d/modules.sh
fi

module load $ATAC_SAMTOOLS_MODULE

for ATAC_TOOL in samtools bedtools; do
    if ! command -v $ATAC_TOOL >/dev/null 2>&1; then
        echo
        echo "    FATAL ERROR: $ATAC_TOOL not on PATH"
        echo
        exit 1
    fi
done
unset ATAC_TOOL
