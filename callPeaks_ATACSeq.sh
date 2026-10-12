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
# From:
#   https://github.com/kundajelab/atac_dnase_pipelines/blob/803e3d92435f52ef31b02fa264fecb5706d550e3/atac.bds#L16
#   Line: 16
#

smooth_window=150

# ./modules/callpeak_macs2_atac.bds:    shiftsize := round( smooth_window.parseReal()/2.0 )

shiftsize=$((smooth_window / 2))

# From paper (Philip, et al, )

pval_thres=0.01

$MACS callpeak \
    -t $TDIR/cleanBED.bed.gz \
    -f BED \
    -n $PREFIX \
    -g $genome \
    -p $pval_thres \
    --nomodel \
    --shift $shiftsize \
    --extsize $smooth_window \
    --call-summits \
    --outdir $ODIR

MACS_ERROR=$?

deactivate
rm -rf $TDIR

exit $MACS_ERROR

