#!/bin/bash

# CMD:
#    bsub -n 1 -q control -o LSF.CTRL/ -J CTRL.ATAC ./PIPE.sh
#

ls ../results/Proj_06125_*/*/align*/*.bam \
    | xargs -n 1 bsub -o LSF.POST/ -J POST2_$$ -R "rusage[mem=24]" -M 25 ../postMapBamProcessing_ATACSeq.sh

bSync POST2_$$
ls *.bed.gz \
    | xargs -n 1 bsub -o LSF.CALLP/ -J CALLP2_$$ -n 3 -R "rusage[mem=16]" -M 17 ../callPeaks_ATACSeq.sh

bSync CALLP2_$$
ls callpeaks/*/*narrow* | fgrep s_N \
    | xargs bsub -o LSF.IDR/ -J IDR -R "rusage[mem=8]" -M 9 ../runIDR.sh

