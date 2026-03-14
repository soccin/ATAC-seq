#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

RESDIR=$1

if [ "$#" != "1" ]; then
    echo
    echo usage: deliveryResults.sh RESDIR
    echo "    RESDIR=/ifs/res/seq/pi/invest/r_###"
    echo
    exit
fi

echo \$RESDIR=$RESDIR $(realpath $RESDIR)

if [ ! -e atacSeq/atlas/macsPeaksMerged.saf ]; then
    echo
    echo ERROR Postprocessing failed
    echo
    exit 1
fi

rsync -avP atacSeq $RESDIR

ATAC_PROJECT_NUM=$($SDIR/extractProjectIDFromPath.py $(realpath $RESDIR) | sed 's/^Proj_//')

cat $SDIR/tmpldeliveryEmail_01.txt \
    | sed 's@ATAC_PROJECT_NUM@'"$ATAC_PROJECT_NUM"'@' \
    | tee DELIVERY_EMAIL_${ATAC_PROJECT_NUM}

DIFFANALYSIS=$(ls | egrep __DiffPeaksEdgeR.*.xlsx)
if [ -e "$DIFFANALYSIS" ]; then

    mkdir -p $RESDIR/atacSeq/diff
    cp *__DiffPeaksEdgeR*.xlsx $RESDIR/atacSeq/diff
    cp *__DiffPeaks*.pdf $RESDIR/atacSeq/diff

    cat $SDIR/tmpldeliveryEmail_02.txt \
        | sed 's@ATAC_PROJECT_NUM@'"$ATAC_PROJECT_NUM"'@' \
        | tee -a DELIVERY_EMAIL_${ATAC_PROJECT_NUM}

fi

cat $SDIR/tmpldeliveryEmail_99.txt \
    | tee -a DELIVERY_EMAIL_${ATAC_PROJECT_NUM}


Rscript ~/Code/BIC/Delivery/Version2j/readme2yaml.R

python3 ~/Code/BIC/Delivery/Version2j/authorization_db/init_impact_project_permissions.py -p project.yaml

