#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ "$#" == "0" ]; then
    echo
    echo "    usage: postMapBamProcessing_ATAC.sh [(-q | --mapq) MAPQ] INPUT_BAM [OUTPUT_BAM]"
    echo
    exit
fi

MAPQ=10

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -q|--mapq)
      MAPQ="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

echo "MAPQ = ${MAPQ}"

IBAM=$1

if [ "$#" == "1" ]; then
    OBAM=${IBAM/.bam/_postProcess.bam}
    OBAM=$(basename $OBAM)
    echo $OBAM
else
    OBAM=$2
fi

TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
#TDIR=_scratch_ATACSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

# f 3 ==> paired, proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
# Also remove reads that have an Insert Size <= 30
#
samtools view -h -q $MAPQ -f 3 -F 1804 $IBAM \
    | awk 'substr($0,1,1)=="@" || sqrt($9*$9)>25' \
    | samtools view -Sb - >$TDIR/step1.bam
picardV2 SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000
picardV2 CollectInsertSizeMetrics LEVEL=null LEVEL=SAMPLE I=$OBAM O=${OBAM/.bam/___INS.txt} H=${OBAM/.bam/___INS.pdf} &

#
# Do Tn5 shift and remove non-standard chromosomes
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | egrep -v "chrUn|_random|GL|NC_|hs37d5|_unplaced" \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shifted.bed.gz}

rm -rf $TDIR

