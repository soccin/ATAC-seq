#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

#
# picard is the ~/bin wrapper (java -Xmx24g). The pipe.sh call site asks
# -M 32G to leave room for the JVM overhead above that heap.
#
if ! command -v picard >/dev/null 2>&1; then
    echo
    echo "    FATAL ERROR: picard not on PATH"
    echo
    exit 1
fi

set -o pipefail

usage() {
    echo
    echo "    usage: postMapBamProcessing_ATACSeq.sh [(-q | --mapq) MAPQ] GENOME INPUT_BAM [OUTPUT_BAM]"
    echo
    echo "    GENOME : b37 | b38 | mm10"
    echo
    exit "${1:-0}"
}

[ "$#" == "0" ] && usage

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

GENOME=$1
shift 1

case $GENOME in
  b37)
  GENOME_BED=$SDIR/lib/genomes/human_b37.genome.bed
  ;;

  b38)
  GENOME_BED=$SDIR/lib/genomes/human_b38.genome.bed
  ;;

  mm10)
  GENOME_BED=$SDIR/lib/genomes/mouse_mm10.genome.bed
  ;;

  *)
  echo -e "\n\nunknown" \$GENOME=$GENOME "\n\n"
  exit 1
  ;;
esac


[ "$#" == "0" ] && usage 1

IBAM=$1

#
# Check that this is a bam file (ie that we did not forget GENOME)
#

samtools quickcheck "$IBAM" || { echo "ERROR: [$IBAM] is not a valid BAM file"; exit 1; }

SID=$(samtools view -H $IBAM | perl -ne '/^@RG.*SM:(\S+)/ && print "$1\n"' | head -1)

if [ "$#" == "1" ]; then

    ODIR=out/$SID
    mkdir -p $ODIR
    OBAM=$ODIR/${SID}_postProcess.bam
    echo \$OBAM=$OBAM

else

    OBAM=$2

fi

#
# Intermediates go to shared scratch (WekaFS, visible from every node).
# /scratch/socci, the JUNO path, is gone.
#
TDIR=${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

# f 3 ==> paired, proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
# Also remove reads that have an Insert Size <= 30
#
samtools view -h -q $MAPQ -f 3 -F 1804 $IBAM \
    | awk 'substr($0,1,1)=="@" || sqrt($9*$9)>30' \
    | samtools view -Sb - >$TDIR/step1.bam \
    || { echo "FATAL ERROR: filter step failed"; exit 1; }

picard SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 \
    || { echo "FATAL ERROR: SortSam failed"; exit 1; }

#
# Runs in the background while the BED is built; waited on below so the
# job does not end, and Slurm does not kill it, before it has finished.
#
picard CollectInsertSizeMetrics LEVEL=null LEVEL=SAMPLE I=$OBAM O=${OBAM/.bam/___INS.txt} H=${OBAM/.bam/___INS.pdf} &
INS_PID=$!

#
# Do Tn5 shift and remove non-standard chromosomes (keep those in genome)
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | bedtools intersect -nonamecheck -a - -b $GENOME_BED \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shifted.bed.gz} \
    || { echo "FATAL ERROR: Tn5 shift / bamtobed failed"; exit 1; }

wait $INS_PID || { echo "FATAL ERROR: CollectInsertSizeMetrics failed"; exit 1; }

rm -rf $TDIR
md5sum ${OBAM/.bam/.shifted.bed.gz} >${OBAM/.bam/.shifted.bed.gz}.md5
