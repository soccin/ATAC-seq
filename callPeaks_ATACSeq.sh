#!/bin/bash

export PYTHONPATH=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/lib/python2.7/site-packages/:$PYTHONPATH
MACS=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/bin/macs2

SDIR="$( cd "$( dirname "$0" )" && pwd )"

IBED=$1

if [ "$#" == "1" ]; then
    PREFIX=${IBED/.bed}
    echo $PREFIX
else
    PREFIX=$2
fi

ODIR=callpeaks/$PREFIX
mkdir -p $ODIR

TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR, $ODIR

export TMPDIR=$TDIR

#
# Clean up BED file to remove non-main chromosomes
#

zcat $IBED \
    | egrep -v "chrUn|_random" \
    | gzip -c - > $TDIR/cleanBED.bed.gz

# MACS2 args
genome=mm

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

exit $MACS_ERROR
