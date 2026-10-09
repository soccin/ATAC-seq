#!/bin/bash
#
# lsfTools.sh -- helper functions for LSF job management
#
# Usage: source lsfTools.sh
#

# bCheck JOBNAME
#
# Call immediately after bSync JOBNAME. Counts jobs in EXIT state
# and aborts the calling script if any failed.
#
bCheck() {
    local jobname=$1
    local n_exit
    n_exit=$(bjobs -noheader -a -J "$jobname" \
             | grep -w "EXIT" | wc -l)
    if [ "$n_exit" -gt 0 ]; then
        echo
        echo "    FATAL: $n_exit job(s) failed for [$jobname]"
        echo
        exit 1
    fi
}
