#!/bin/bash
#
# genomeInfo.sh TAG FIELD
#
# Print FIELD of genome TAG from lib/genomes/GENOMES.tsv. A supported
# genome is a row of that table; FIELD is one of its column names
# (TAG, SPECIES, MACS_GSIZE, TXDB, ORGDB, ADD_CHR, MITO, MSIGDB_DB,
# MSIGDB_SETS). An unknown tag
# or field is an error: the message lists what the table has, and the
# exit code is 1.
#
# genomeInfo.sh TAG TAG prints TAG if it is supported, so it doubles as
# the supported-genome check.
#
# The genome file of TAG, lib/genomes/TAG.genome, lists the contigs
# analyzed (autosomes, X and Y) with their lengths; it is the region
# filter and the chromosome sizes for every stage.
#

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"
TABLE=$SDIR/../lib/genomes/GENOMES.tsv

if [ "$#" != "2" ]; then
    echo "usage: genomeInfo.sh TAG FIELD" >&2
    exit 1
fi

TAG=$1
FIELD=$2

awk -F'\t' -v tag="$TAG" -v field="$FIELD" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
            fields = fields (i > 1 ? " " : "") $i
        }
        next
    }
    {
        tags = tags (tags != "" ? " " : "") $1
        if ($1 == tag) {
            row = $0
        }
    }
    END {
        if (!(field in col)) {
            printf "genomeInfo.sh: unknown field [%s]; fields: %s\n", field, fields >"/dev/stderr"
            exit 1
        }
        if (row == "") {
            printf "genomeInfo.sh: unsupported genome [%s]; supported: %s\n", tag, tags >"/dev/stderr"
            exit 1
        }
        split(row, f, "\t")
        print f[col[field]]
    }' "$TABLE"
