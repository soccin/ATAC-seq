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
    --outdir $ODIR

MACS_ERROR=$?

deactivate

if [ "$MACS_ERROR" != "0" ]; then
    echo "FATAL ERROR: macs2 callpeak failed (rc $MACS_ERROR)"
    rm -rf $TDIR
    exit 1
fi

#
# ENCODE's steps after MACS2 (encode_task_macs2_atac.py), then its
# blacklist filter (encode_lib_blacklist_filter.py, run on that result
# by encode_task_post_call_peak_atac.py):
#
#   1. sort by -log10 p (column 8), highest first
#   2. clamp a negative start or end to 0; a missing summit (column 10
#      -1) becomes the peak center
#   3. keep the top 300,000 peaks
#   4. clip to the contig ends: UCSC bedClip -truncate, the v369 build
#      ENCODE installs (scripts/install_conda_env.sh)
#   5. drop each peak that overlaps a blacklisted region and cap the
#      score (column 5) at 1000
#
# ENCODE also renames the peaks Peak_1, Peak_2, ...; here they keep
# MACS2's names, so each still matches its summit and its .xls row. In
# step 2 ENCODE writes the peak center as a contig position; column 10
# is an offset from the start, so it is written as one here. MACS2 with
# --call-summits does not write -1, so neither form has been used.
#
# MACS2's own files are left as they are. The result, and the summits
# of the peaks kept (matched by peak name), are written beside them as
# *.blacklistFiltered.*. The atlas, and so the counts, use only the
# blacklistFiltered peaks.
#

MAX_PEAKS=300000

RAW_PEAKS=$ODIR/${PREFIX}_peaks.narrowPeak
RAW_SUMMITS=$ODIR/${PREFIX}_summits.bed
BL_PEAKS=$ODIR/${PREFIX}_peaks.blacklistFiltered.narrowPeak
BL_SUMMITS=$ODIR/${PREFIX}_summits.blacklistFiltered.bed

for F in $RAW_PEAKS $RAW_SUMMITS; do
    if [ ! -e "$F" ]; then
        echo "FATAL ERROR: macs2 did not write $F"
        rm -rf $TDIR
        exit 1
    fi
done

LC_ALL=C sort -k8,8gr -T "$TMPDIR" $RAW_PEAKS \
  | awk 'BEGIN {OFS = "\t"}
         {
             if ($2 < 0) $2 = 0
             if ($3 < 0) $3 = 0
             if ($10 == -1) $10 = int(($3 - $2 + 1) / 2)
             print
         }' \
  > $TDIR/sorted.narrowPeak \
  || { echo "FATAL ERROR: sort of $RAW_PEAKS failed"; rm -rf $TDIR; exit 1; }

#
# head on its own: in a pipe, sort and awk would be stopped by SIGPIPE
# when head exits, and pipefail would report that as a failure.
#
head -n $MAX_PEAKS $TDIR/sorted.narrowPeak > $TDIR/top.narrowPeak \
  || { echo "FATAL ERROR: cap of $RAW_PEAKS failed"; rm -rf $TDIR; exit 1; }

$SDIR/bin/bedClip $TDIR/top.narrowPeak $GENOME_FILE $TDIR/clipped.narrowPeak \
    -truncate -verbose=2 \
  || { echo "FATAL ERROR: bedClip of $RAW_PEAKS failed"; rm -rf $TDIR; exit 1; }

#
# bedClip -truncate turns a peak that lies wholly past the contig end
# into start > end, and exits 0. MACS2 cannot call one (a peak holds
# reads on the contig); stop if it ever happens.
#
NBAD=$(awk '$2 >= $3' $TDIR/clipped.narrowPeak | wc -l)
if [ "$NBAD" != "0" ]; then
    echo "FATAL ERROR: $NBAD peaks of $RAW_PEAKS are past their contig end"
    rm -rf $TDIR
    exit 1
fi

bedtools intersect -nonamecheck -v -a $TDIR/clipped.narrowPeak -b $BLACKLIST \
  | awk 'BEGIN {OFS = "\t"} {if ($5 > 1000) $5 = 1000; print}' \
  > $BL_PEAKS \
  || { echo "FATAL ERROR: blacklist filter of $RAW_PEAKS failed"; rm -rf $TDIR; exit 1; }

awk 'NR == FNR {keep[$4] = 1; next} $4 in keep' $BL_PEAKS $RAW_SUMMITS \
  > $BL_SUMMITS \
  || { echo "FATAL ERROR: summit filter failed"; rm -rf $TDIR; exit 1; }

echo "peaks of $PREFIX: MACS2 $(wc -l <$RAW_PEAKS)," \
     "top $(wc -l <$TDIR/top.narrowPeak) (cap $MAX_PEAKS)," \
     "kept $(wc -l <$BL_PEAKS) after blacklist $BLACKLIST"

rm -rf $TDIR

exit 0

