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

#
# Temp files never go to /tmp. IRIS compute nodes set TMPDIR=/tmp, a
# 137G volume shared by every job on the node; sort spill, the R session
# tempdir and Python tempfiles all land there by default. Point TMPDIR at
# the node-local /localscratch (2.8T) instead. bin/runStage.sh has
# already set a per-job directory under it for stage jobs; keep that.
# This also covers stage scripts run by hand.
#

ATAC_LOCAL_TMP=${ATAC_LOCAL_TMP:-/localscratch/$USER}

if [[ "$TMPDIR" != "$ATAC_LOCAL_TMP"/* ]]; then
    export TMPDIR=$ATAC_LOCAL_TMP
fi

if ! mkdir -p "$TMPDIR"; then
    echo
    echo "    FATAL ERROR: cannot create TMPDIR=[$TMPDIR]"
    echo
    exit 1
fi
