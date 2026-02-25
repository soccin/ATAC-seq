#!/bin/bash

module load samtools

if [ "$#" != "1" ]; then
    echo usage getGenomeBuild.sh BAM
    exit
fi

SAMTOOLS=$(which samtools)
if [ $SAMTOOLS == "" ]; then
    echo samtools not in current path
    exit -1
fi

GENOME_MD5=$($SAMTOOLS view -H $1 | egrep "^@SQ" | cut -f-3 | sort  | md5sum - | awk '{print $1}')

case $GENOME_MD5 in
    b879c678e7fd80718cf20d10c6b846e4)
    # b37 gatk /ifs/depot/assemblies/H.sapiens/b37/b37.dict
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

    7f8c5323ff7e0ff6b5d40efe53eaf050)
    # BIC Xeno-graph genome
    echo "b37+mm10"
    ;;

    d660fd17a979374182d3ba8b6d76cac0)
    # UCSC mm10 /ifs/depot/assemblies/M.musculus/mm10/mm10.dict
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

