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
    echo "    GENOME : a tag in lib/genomes/GENOMES.tsv (b37 | hg38 | mm10)"
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

#
# lib/genomes/<GENOME>.genome lists the contigs analyzed: autosomes, X
# and Y. Mitochondrial, unplaced, random, decoy and viral contigs are
# left out of the shifted BED, so out of peaks, bigWigs and the atlas.
#
MITO=$($SDIR/bin/genomeInfo.sh "$GENOME" MITO) || exit 1
GENOME_FILE=$SDIR/lib/genomes/$GENOME.genome
if [ ! -s "$GENOME_FILE" ]; then
    echo "FATAL ERROR: no genome file $GENOME_FILE"
    exit 1
fi


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
#
# ENCODE's filter (encode_task_filter.py): MAPQ and flags, then
# samtools fixmate -r on the name-sorted reads, which clears the proper
# pair flag of a read whose mate was removed, and the flag filter again,
# which removes it. Until v2.0.0 such orphans were kept.
#
# Two filters ENCODE does not have. Both remove reads that are not from
# the reference. BWA-MEM is a local aligner: it aligns the best 19-35
# bp seed of such a read, soft-clips the rest and gives the seed a MAPQ
# of its own. ENCODE aligns end to end with bowtie2, which leaves these
# reads unaligned. In Proj_18143 they were a tandem repeat, the same in
# every sample, up to 55% of a sample's reads, piled up by the hundreds
# on single bases, and they made up 38-75% of the counts of ten of its
# 15 samples. MarkDuplicates cannot see them: the clip length varies,
# which moves the unclipped 5' position it compares.
#
# KEEP_5P, on the first view: the 5' end of a read must align, as a
# Tn5 read starts at the cut site. The 5' end is the start of a forward
# read and the end of a reverse read (0x10; CIGAR is in reference
# order). A 3' clip is kept: it is mostly the Nextera adapter, which
# PEMapper does not trim. It is applied before the name sort so that
# fixmate -r and the second -f 3 remove the mate of a dropped read.
#
# The insert filter, in this pipeline since 2022: pairs with
# |TLEN| <= 30 are removed (fixmate gives both reads the same |TLEN|).
# After KEEP_5P it removes the few seed pairs left whose 5' ends both
# align, e.g. 22M79S with its mate 79S22M.
#
KEEP_5P='(flag.reverse && !(cigar =~ "S$")) || (!flag.reverse && !(cigar =~ "^[0-9]+S"))'

samtools view -u -q $MAPQ -f 3 -F 1804 -e "$KEEP_5P" $IBAM \
    | samtools sort -n -m 4G -T $TMPDIR/nsort -u - \
    | samtools fixmate -r - - \
    | samtools view -b -f 3 -F 1804 -e 'tlen > 30 || tlen < -30' \
        -o $TDIR/step1.bam - \
    || { echo "FATAL ERROR: filter step failed"; exit 1; }

picard SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 \
    || { echo "FATAL ERROR: SortSam failed"; exit 1; }

#
# Runs in the background while the BED is built; waited on below so the
# job does not end, and Slurm does not kill it, before it has finished.
#
# W=1000 fixes the histogram at 0-1000 bp, as ENCODE's fraglen_stat_pe
# does and plotINSStats.R plots. Without W, picard cuts the histogram
# (and the mean and SD) at median + 10 x MAD; in a sample with many
# short fragments that falls inside the nucleosome peaks.
#
picard CollectInsertSizeMetrics LEVEL=null LEVEL=SAMPLE W=1000 \
    I=$OBAM O=${OBAM/.bam/___INS.txt} H=${OBAM/.bam/___INS.pdf} &
INS_PID=$!

#
# Do Tn5 shift and keep only the contigs of the genome file. On the way,
# count the reads on the mitochondrial contig ($MITO) out of all reads
# in the filtered BAM, before they are dropped, for the ___MITO.txt
# metric.
#

awk 'BEGIN {OFS = "\t"} {print $1, 0, $2}' $GENOME_FILE >$TDIR/genome.bed \
    || { echo "FATAL ERROR: cannot write $TDIR/genome.bed"; exit 1; }

MITO_TXT=${OBAM/.bam/___MITO.txt}

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | awk -F'\t' -v mito="$MITO" -v sid="$SID" -v out="$MITO_TXT" '
        { n++; if ($1 == mito) m++; print }
        END {
            printf "SAMPLE\tMITO_CONTIG\tREADS\tMITO_READS\tFRAC_MITO\n" >out
            printf "%s\t%s\t%d\t%d\t%.4f\n", sid, mito, n, m, (n ? m / n : 0) >out
        }' \
    | bedtools intersect -nonamecheck -a - -b $TDIR/genome.bed \
    | awk -F'\t' \
        'BEGIN {OFS = FS} { if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
    | gzip -nc >${OBAM/.bam/.shifted.bed.gz} \
    || { echo "FATAL ERROR: Tn5 shift / bamtobed failed"; exit 1; }

wait $INS_PID || { echo "FATAL ERROR: CollectInsertSizeMetrics failed"; exit 1; }

rm -rf $TDIR
md5sum ${OBAM/.bam/.shifted.bed.gz} >${OBAM/.bam/.shifted.bed.gz}.md5
