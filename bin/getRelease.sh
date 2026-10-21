#!/bin/bash

#
# getRelease.sh -- print the release of this checkout, the
# "- **Version:** vX.Y.Z" line of VERSION.md
#
# The release is the client-facing version: it is set by hand in
# VERSION.md when a release is made, never derived from git. pipe.sh
# records it as RELEASE= in 00.RUNSTATUS.txt; deliverResults.sh puts
# it in the email and the doc links.
#

set -euo pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"
VFILE=$(dirname "$SDIR")/VERSION.md

RELEASE=$(sed -n 's/^- \*\*Version:\*\* *//p' "$VFILE" 2>/dev/null || true)
if [[ ! "$RELEASE" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "getRelease.sh: no '- **Version:** vX.Y.Z' line in $VFILE" >&2
    exit 1
fi

echo "$RELEASE"
