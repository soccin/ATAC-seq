#!/bin/bash


SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh
. $SDIR/venv/bin/activate
MACS=macs2

set -o pipefail

GENOMEBUILD=$1
IBED=$2

genome=$($SDIR/bin/genomeInfo.sh "$GENOMEBUILD" MACS_GSIZE) || exit 1
GENOME_FILE=$SDIR/lib/genomes/$GENOMEBUILD.genome
if [ ! -s "$GENOME_FILE" ]; then
    echo "FATAL ERROR: no genome file $GENOME_FILE"
    exit 1
fi
BLACKLIST=$SDIR/lib/genomes/$GENOMEBUILD.blacklist.bed
if [ ! -s "$BLACKLIST" ]; then
    echo "FATAL ERROR: no blacklist $BLACKLIST"
    exit 1
fi

echo "md5sum ${IBED}"
echo $(md5sum ${IBED})

checkMD5=$(md5sum -c ${IBED}.md5 | awk '{print $2}'); 
if [ "$checkMD5" != "OK" ]; then 
    echo
    echo md5check fail
    echo sleeping for 300 secs 
    sleep 300
    checkMD5=$(md5sum -c ${IBED}.md5 | awk '{print $2}'); 
    if [ "$checkMD5" != "OK" ]; then 
        echo 
        echo Second md5check fail
        echo FATAL ERROR Stopping
        echo
        exit 1
    fi
fi


if [ "$#" == "2" ]; then
    PREFIX=$(basename ${IBED/.bed.*})
else
    PREFIX=$3
fi
echo \$PREFIX=$PREFIX

ODIR=callpeaks/$PREFIX
mkdir -p $ODIR

TDIR=${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR, $ODIR

export TMPDIR=$TDIR

#
# Clean up BED file to keep only the contigs of the genome file
#

awk 'BEGIN {OFS = "\t"} {print $1, 0, $2}' $GENOME_FILE >$TDIR/genome.bed \
  || { echo "FATAL ERROR: cannot write $TDIR/genome.bed"; rm -rf $TDIR; exit 1; }

zcat $IBED \
  | bedtools intersect -nonamecheck -a - -b $TDIR/genome.bed \
  | gzip -c - > $TDIR/cleanBED.bed.gz \
  || { echo "FATAL ERROR: chromosome filter of $IBED failed"; rm -rf $TDIR; exit 1; }

#
# MACS2 settings of the ENCODE ATAC-seq pipeline,
# src/encode_task_macs2_atac.py:
#
#   shiftsize = -int(round(float(smooth_win)/2.0))
#   macs2 callpeak -t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh}
#     --shift {shiftsize} --extsize {extsize}
#     --nomodel -B --SPMR --keep-dup all --call-summits
#
# -B --SPMR only write signal bedGraphs, which this pipeline does not use
# (its bigWigs come from makeBigWigFromBEDZ.sh), so they are left out.
#

smooth_window=150

#
# Negative: move each cut site (the 5' end of the shifted read) back half
# the window, so --extsize centers a smooth_window-wide region on it.
# Until v2.0.0 this was +75, which put the region 75-225 bp downstream of
# the cut site.
#
shiftsize=$(( -smooth_window / 2 ))

pval_thres=0.01

#
# --keep-dup all: duplicates are already removed in POST (-F 1804), so
# reads with the same cut site are distinct fragments. Until v2.0.0 the
# MACS2 default (1) kept one of them.
#
# MACS2 writes to scratch; only the blacklist-filtered peaks below go to
# callpeaks/.
#
mkdir -p $TDIR/macs

$MACS callpeak \
    -t $TDIR/cleanBED.bed.gz \
    -f BED \
    -n $PREFIX \
    -g $genome \
    -p $pval_thres \
    --nomodel \
    --shift $shiftsize \
    --extsize $smooth_window \
    --keep-dup all \
    --call-summits \
    --outdir $TDIR/macs

MACS_ERROR=$?

deactivate

if [ "$MACS_ERROR" != "0" ]; then
    echo "FATAL ERROR: macs2 callpeak failed (rc $MACS_ERROR)"
    rm -rf $TDIR
    exit 1
fi

#
# Blacklist filter, as ENCODE's encode_lib_blacklist_filter.py does for
# every peak set: drop each peak that overlaps a blacklisted region and
# cap the score (column 5) at 1000. Only the filtered peaks are written
# to callpeaks/, so the atlas, the counts and the delivered narrowPeak
# all leave them out. The summits are those of the peaks kept (matched
# by peak name). The .xls is MACS2's own table, before the filter.
#

RAW=$TDIR/macs/${PREFIX}_peaks.narrowPeak

bedtools intersect -nonamecheck -v -a $RAW -b $BLACKLIST \
  | awk 'BEGIN {OFS = "\t"} {if ($5 > 1000) $5 = 1000; print}' \
  > $ODIR/${PREFIX}_peaks.narrowPeak \
  || { echo "FATAL ERROR: blacklist filter of $RAW failed"; rm -rf $TDIR; exit 1; }

awk 'NR == FNR {keep[$4] = 1; next} $4 in keep' \
    $ODIR/${PREFIX}_peaks.narrowPeak $TDIR/macs/${PREFIX}_summits.bed \
  > $ODIR/${PREFIX}_summits.bed \
  || { echo "FATAL ERROR: summit filter failed"; rm -rf $TDIR; exit 1; }

cp $TDIR/macs/${PREFIX}_peaks.xls $ODIR \
  || { echo "FATAL ERROR: cannot copy ${PREFIX}_peaks.xls"; rm -rf $TDIR; exit 1; }

echo "blacklist $BLACKLIST: kept $(wc -l <$ODIR/${PREFIX}_peaks.narrowPeak)" \
     "of $(wc -l <$RAW) peaks"

rm -rf $TDIR

exit 0

