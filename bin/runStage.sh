#!/bin/bash
#
# runStage.sh CMD [ARG ...]
#
# Run one stage command and append an exit trailer to the job log.
# pipe.sh submits every stage job through this script.
#
# Slurm writes nothing of its own into a job's -o file, so a log alone
# cannot show whether the job worked. The trailer
#
#   #ATAC_EXIT=<rc>
#
# is the last line the job writes. A missing trailer means the job never
# reached the end (killed for memory or time, node failure); then sacct
# is the only record. bin/checkRun.sh reads the trailer as a fallback when
# sacct no longer has the job.
#
# The trailer goes to stderr. The MergePeaks stage redirects the command's
# stdout into macsPeaksMerged.saf; stderr still lands in the job log.
#

echo "#ATAC_HOST=$(hostname) START=$(date '+%Y-%m-%d %H:%M:%S')" >&2
echo "#ATAC_CMD=$*" >&2

"$@"
RC=$?

echo "#ATAC_EXIT=$RC END=$(date '+%Y-%m-%d %H:%M:%S')" >&2
exit $RC
