#!/bin/bash


SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/venv/bin/activate
MACS=macs2

GENOMEBUILD=$1
IBED=$2

if [ "$#" == "2" ]; then
    PREFIX=${IBED/.bed}
    echo $PREFIX
else
    PREFIX=$3
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
    | egrep -v "chrUn|_random|_unplaced" \
    | gzip -c - > $TDIR/cleanBED.bed.gz

# MACS2 args

case $GENOMEBUILD in

    b37)
    genome=hs
    ;;

    mm10)
    genome=mm
    ;;

    sCer+sMik_IFO1815)
    genome=23606800
    ;;

    *)
    echo
    echo "    Unknown GENOMEBUILD [$GENOMEBUILD]"
    echo
    exit 1
    ;;

esac

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

exit $MACS_ERROR

