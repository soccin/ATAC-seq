#!/bin/bash
#
# getGenomeBuildBAM.sh BAM
#
# Print the genome build tag of BAM, identified by the md5 of its sorted
# @SQ SN/LN lines. Prints "unknown <md5>" for a header not in the table.
# Exits 1 if the header cannot be read.
#
# The table names more builds than the pipeline supports, so that an
# unsupported BAM gets a named error rather than "unknown". The
# supported builds are the rows of lib/genomes/GENOMES.tsv.
#
# Keep this file identical in ChIP-seq and ATAC-seq.
#

set -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/loadTools.sh

if [ "$#" != "1" ]; then
    echo usage getGenomeBuildBAM.sh BAM
    exit 1
fi

HEADER=$(samtools view -H "$1") \
    || { echo "FATAL ERROR: cannot read the header of [$1]" >&2; exit 1; }

GENOME_MD5=$(echo "$HEADER" | egrep "^@SQ" | cut -f-3 | sort | md5sum - | awk '{print $1}') \
    || { echo "FATAL ERROR: cannot compute the @SQ md5 of [$1]" >&2; exit 1; }

case $GENOME_MD5 in
    b879c678e7fd80718cf20d10c6b846e4)
    # b37 gatk; IRIS H.sapiens/b37 and b37_gatk (PEMapper human_b37)
    echo "b37"
    ;;

    117fce86b797081e0af6d69cbd94dcde)
    # b37 version used by DMP pipeline
    echo "b37_dmp"
    ;;

    5b4e380a6b4fc3494cfc66c917d41b37)
    # UCSC hg19 /ifs/depot/assemblies/H.sapiens/hg19/hg19.dict
    echo "hg19"
    ;;

    3d72c6961689390556ed2d5a33e66e17)
    # Main chromosomes only (used by cfDNA collaboration)
    echo "hg19-mainOnly"
    ;;

    933b376d936c265fc6b44c8bd19fc66d)
    # TCGA BAMs UR:ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz
    # AS:GRCh37-lite (b37-ish)
    echo "GRCh37-lite"
    ;;

    6a02a9445433ae0b01febe3b4f39d5ec)
    # GRCh38: IRIS H.sapiens/GRCh38.d1.vd1, the GDC reference (chr1, chr2, ...)
    # (PEMapper human_GRCh38). Tagged hg38, the name the Broad and UCSC use;
    # its main contigs equal UCSC hg38's in name and length.
    echo "hg38"
    ;;

    7f8c5323ff7e0ff6b5d40efe53eaf050)
    # BIC Xeno-graph genome
    echo "b37+mm10"
    ;;

    863ad3aa4e086b3224adce2aa4f48b01)
    # mouse_hybrid (MM10+B37)
    echo "mm10+b37"
    ;;

    a10cc06a54dcc732a20cee167d1f1afe)
    # hybrids/H.sapiens_M.musculus/hg38_mm39 (PEMapper human_hg38+mm39)
    echo "hg38+mm39"
    ;;

    d660fd17a979374182d3ba8b6d76cac0)
    # UCSC mm10; IRIS M.musculus/mm10 (PEMapper mouse_mm10)
    echo "mm10"
    ;;

    34839afd79d8b772037ed7a0e0a4f9c3)
    # UCSC mm10
    echo "mm10_hBRAF_V600E"
    ;;

    f9cd233a3d5c9540eece434c65f84f1c)
    # mm9 Full
    echo "mm9Full"
    ;;

    0835b244e0adb20253e5fa1c4ee58ec4)
    # mouse_GRCm38
    echo "GRC_m38"
    ;;

    8a300152df87118834c4268ab4b713aa)
    # Yeast hybrid sCer+sMik_IFO1815
    echo "sCer+sMik_IFO1815"
    ;;

    *)
    echo "unknown" $GENOME_MD5
    ;;
esac
