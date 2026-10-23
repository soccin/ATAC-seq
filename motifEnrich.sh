#!/bin/bash

#
# motifEnrich.sh -- motif enrichment of the differential peaks, as a job
#
# usage, from the analysis directory, after R/diffAnalysisPairwise.R:
#   bsub -n 8 -M 128G -W 6:00:00 -J motifs -o SLURM.MOTIF/ \
#       /path/to/ATAC-seq/motifEnrich.sh [-w WIDTH] [-m MINPEAKS] [RUNTAG]
#
# Runs R/motifEnrich.R with the tools of bin/loadTools.sh, which also
# keeps R's temp files out of /tmp. The motif scan uses every core of
# the job (bsub -n sets SLURM_CPUS_PER_TASK, the default of -c).
# Outside a job, or with -h, it prints the bsub line and how to choose
# the options and the job size (usage below). README.md, "Motif
# enrichment", has the same notes; keep the two the same.
# R/motifEnrich.R --help describes the method and the output.
#

set -euo pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd -P )"

usage() {
    cat <<EOF

usage: motifEnrich.sh [-w WIDTH] [-m MINPEAKS] [RUNTAG]

Run as a job from the analysis directory, after diffAnalysisPairwise.R,
with the same RUNTAG (default: none):

  bsub -n 8 -M 128G -W 6:00:00 -J motifs -o SLURM.MOTIF/ \\
      $SDIR/motifEnrich.sh [-w WIDTH] [-m MINPEAKS] [RUNTAG]

Writes <PROJ>_<RUNTAG>_MotifEnrichment.xlsx, .pdf and .rds to
atacSeq/diff; the log is SLURM.MOTIF/<jobid>.out.

-w WIDTH (default 501; odd, at least 51)
    The bp scanned for motifs around each peak's strongest summit, for
    the up, down and background peaks alike. 501 (summit +/- 250) is
    ArchR's peak width and in chromVAR's 250-500 range for ATAC, which
    it calls "fairly robust to the exact choice". HOMER's default,
    from TF ChIP-seq, is 200. Keep 501 unless you have a reason: 201
    (summit +/- 100) keeps closer to the summit; a wider window adds
    more flanking sequence. The scan time grows with the width.

-m MINPEAKS (default 50)
    An up or down set with fewer peaks is not tested; the Summary sheet
    says so. 50 is a judgment, not a standard: each set is tested for
    879 motifs, with BH. In Proj_18143 ED1 vs EV (33 up, 17 down) was
    not tested; of the sets tested, 244 peaks gave no motif with
    padj < 0.05 and 441 gave 92. Lower it to try a weak comparison and
    read its hits with care; raising it only skips more sets.

Job size
    -n  cores for the motif scan, about 26 sequences per second per
        core. Each comparison with a set tested scans all its tested
        peaks, so the time is about
          (sum of tested peaks over those comparisons) / (26 * n) s.
        A comparison's tested peaks are the "N peaks" on its MA and
        volcano plots in <PROJ>_<RUNTAG>_DiffPeaks_V5.pdf. Proj_18143:
        6 comparisons scanned, 64,000 to 240,000 peaks each, 1.04
        million in all, 80 min with -n 8.
    -M  Proj_18143 reached 64G with -n 8; 128G leaves room. With more
        cores give more; how much has not been measured.
    -W  set from the time above, with room; 6:00:00 for a run like
        Proj_18143.

The same notes: $SDIR/README.md,
"Motif enrichment", "Choosing the options and the job size".
The method and output: Rscript $SDIR/R/motifEnrich.R --help

EOF
}

for arg in "$@"; do
    if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
        usage
        exit 0
    fi
done

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo
    echo "    motifEnrich.sh must run as a job; submit it with bsub"
    usage
    exit 1
fi

. "$SDIR/bin/loadTools.sh"

echo "JOB=$SLURM_JOB_ID HOST=$(hostname) CPUS=${SLURM_CPUS_PER_TASK:-1}"
echo "TMPDIR=$TMPDIR"
echo "PWD=$(pwd -P)"
date
Rscript "$SDIR/R/motifEnrich.R" "$@"
date
