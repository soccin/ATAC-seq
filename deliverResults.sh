#!/bin/bash

#
# deliverResults.sh -- deliver a finished run to the BIC delivery area
#
# usage:
#   deliverResults.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]
#   deliverResults.sh [-n] -d|--default
#
# Run from the analysis directory, after the run has COMPLETED
# (bin/checkRun.sh exits 0). It
#
#   - writes project.yaml for the bicdelivery permissions database,
#     from the run and the first README.txt found going up from here
#     to Proj_<no> (bin/readme2yaml.R);
#   - copies atacSeq/ to DELIVERY/r_NNN/, with the client docs
#     (docs/output.pdf, docs/methods.pdf of this checkout) in
#     atacSeq/docs and the differential analysis files of the project
#     (Proj_<no>_*_DiffPeaksEdgeR_V4.xlsx, Proj_<no>_*_DiffPeaks_V4.pdf),
#     if any, in atacSeq/diff;
#   - writes the delivery email to deliveryEmail_<no>_r_NNN_<yymmdd>.txt;
#   - registers the run with the permissions database.
#
#   -d, --default   deliver to DROOT/<pi>/<invest>/Proj_<no>, taken
#                   from the path .../Users/<Pi>/<Invest>/Proj_<no>/...
#                   of the analysis directory (lowercased but for the
#                   project)
#   -n, --dry-run   print the commands that write to the delivery area
#                   or the database instead of running them. The email
#                   and project.yaml are still written here.
#
# Without a trailing r_NNN the run folder is the one after the highest
# r_NNN in the delivery folder, else r_001. An existing r_NNN can be
# given to add files to it, e.g. a differential analysis made later;
# nothing in it is deleted.
#
# DROOT is /data1/core002/res/bic/results. Every write to the delivery
# area runs as the service account svc_core001_bot01 (dzdo).
#

set -euo pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

BOT_USER=svc_core001_bot01
DROOT=/data1/core002/res/bic/results
BIC_DELIVERY=$HOME/Code/BIC/Delivery/Version2j

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

WARNINGS=()
warn() {
    echo "    WARNING: $*"
    WARNINGS+=("$*")
}

#
# The physical path, as R's getwd() gives it, so that PROJNO is the
# project the R reports were named by.
#
CWD=$(pwd -P)

PROJ=""
IFS=/ read -r -a PARTS <<< "$CWD"
for p in "${PARTS[@]}"; do
    if [[ "$p" == Proj_* ]]; then
        PROJ=$p
    fi
done

#
# Default delivery path from .../Users/Aa/BBB/Proj_nnnnn/...
# -> aaa/bbb/Proj_nnnnn; the third component must be the project.
#
DPATH=""
if [[ "$CWD" == *"/Users/"* && -n "$PROJ" ]]; then
    rest="${CWD#*/Users/}"
    IFS=/ read -r u1 u2 proj _ <<< "$rest"
    if [[ -n "$u1" && -n "$u2" && "$proj" == "$PROJ" ]]; then
        DPATH="${u1,,}/${u2,,}/$proj"
    fi
fi

usage() {
    echo
    echo "usage: deliverResults.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]"
    echo "       deliverResults.sh [-n] -d|--default"
    echo
    echo "    -n, --dry-run   print the delivery and database commands, do not run them"
    echo
    echo "    default: ${DPATH:+$DROOT/}${DPATH:-(none)}"
    echo
}

ODIR=""
USE_DEFAULT=0
DRYRUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--default)
            USE_DEFAULT=1
            ;;
        -n|--dry-run)
            DRYRUN=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [[ -n "$ODIR" ]]; then
                echo "Only one delivery folder may be given"
                usage
                exit 1
            fi
            ODIR="$1"
            ;;
    esac
    shift
done

if [[ "$USE_DEFAULT" -eq 1 && -n "$ODIR" ]]; then
    echo "Give a delivery folder or -d, not both"
    usage
    exit 1
fi

if [[ "$USE_DEFAULT" -eq 1 ]]; then
    if [[ -z "$DPATH" ]]; then
        fatal "no default delivery path for $CWD (not .../Users/<Pi>/<Invest>/Proj_<no>/...)"
    fi
    ODIR="$DROOT/$DPATH"
fi

if [[ -z "$ODIR" ]]; then
    usage
    exit 1
fi

asbot() {
    if [[ "$DRYRUN" -eq 1 ]]; then
        echo "DRYRUN: dzdo -u $BOT_USER ${*@Q}"
    else
        dzdo -u "$BOT_USER" "$@"
    fi
}

##
# The run

STATUSF=00.RUNSTATUS.txt
if [ ! -e "$STATUSF" ]; then
    fatal "no $STATUSF; run deliverResults.sh from the analysis directory"
fi

if ! "$SDIR/bin/checkRun.sh"; then
    fatal "the run has not COMPLETED (bin/checkRun.sh above); nothing delivered"
fi
echo

getKey() {
    sed -n "s/^$1=//p" "$STATUSF" | tail -1
}

GENOME=$(getKey GENOME)
MAPQ=$(getKey MAPQ)
RUN_VERSION=$(getKey VERSION)
for key in GENOME MAPQ RUN_VERSION; do
    if [[ -z "${!key}" ]]; then
        fatal "no ${key#RUN_}= in $STATUSF"
    fi
done

if [[ -z "$PROJ" ]]; then
    fatal "no Proj_* directory in $CWD; cannot name the project"
fi
PROJNO=${PROJ#Proj_}

#
# The email names the release the run was made with: VERSION is
# git describe --long, TAG-N-gSHA. Anything but a vX.Y.Z tag with N=0
# needs the version and the doc links in the email fixed by hand.
#
if [[ "$RUN_VERSION" =~ ^(.+)-([0-9]+)-g[0-9a-f]+$ ]]; then
    VERSION=${BASH_REMATCH[1]}
    NAHEAD=${BASH_REMATCH[2]}
else
    VERSION=$RUN_VERSION
    NAHEAD=""
fi
if [[ ! "$VERSION" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ || "$NAHEAD" != "0" ]]; then
    warn "run version $RUN_VERSION is not a release; fix the version and the doc links in the email"
fi

#
# pipe.sh writes VERSION with the same git describe flags.
#
if ! SDIR_VERSION=$(git -C "$SDIR" describe --always --long); then
    warn "cannot get the version of $SDIR; the docs delivered may not match run $RUN_VERSION"
elif [[ "$SDIR_VERSION" != "$RUN_VERSION" ]]; then
    warn "the docs delivered are from $SDIR_VERSION ($SDIR), the run is $RUN_VERSION"
fi

DOCS=("$SDIR/docs/output.pdf" "$SDIR/docs/methods.pdf")
for f in "${DOCS[@]}"; do
    if [ ! -s "$f" ]; then
        fatal "missing $f; render it with docs/MAKE_DOCS.sh"
    fi
done

#
# Differential analysis: every RUNTAG of this project, untagged
# (Proj_<no>__...) or tagged (Proj_<no>_<tag>_...).
#
shopt -s nullglob
DIFF_XLSX=("${PROJ}"_*_DiffPeaksEdgeR_V4.xlsx)
DIFF_PDF=("${PROJ}"_*_DiffPeaks_V4.pdf)
shopt -u nullglob
DIFF_FILES=("${DIFF_XLSX[@]}" "${DIFF_PDF[@]}")

if [[ ${#DIFF_XLSX[@]} -ne ${#DIFF_PDF[@]} ]]; then
    warn "${#DIFF_XLSX[@]} DiffPeaksEdgeR_V4.xlsx but ${#DIFF_PDF[@]} DiffPeaks_V4.pdf files"
fi
for f in "${DIFF_FILES[@]}"; do
    if [[ "$f" != "${PROJ}__DiffPeaksEdgeR_V4.xlsx" && "$f" != "${PROJ}__DiffPeaks_V4.pdf" ]]; then
        warn "diff file $f has a run tag; the email names the untagged files, edit it"
    fi
done

##
# Run folder: a trailing r_NNN as given, else the next one in ODIR

ODIR=$(realpath -m "$ODIR")
if [[ "$ODIR" =~ /(r_[0-9]+)$ ]]; then
    CURRDIR=${BASH_REMATCH[1]}
    ODIR=${ODIR%/*}
else
    CURRDIR=""
fi

#
# The deepest existing folder of ODIR must be readable by us; else
# the run folders in it could exist unseen, and r_001 would be picked
# or files added to a run folder without the note below.
#
seen=$ODIR
while [[ -n "$seen" && ! -e "$seen" ]]; do
    seen=${seen%/*}
done
seen=${seen:-/}
if [[ ! -r "$seen" || ! -x "$seen" ]]; then
    fatal "cannot read $seen, so cannot tell what is in $ODIR"
fi

if [[ -z "$CURRDIR" ]]; then
    last=0
    for d in "$ODIR"/r_[0-9]*; do
        n=${d##*/r_}
        if [[ -d "$d" && "$n" =~ ^[0-9]+$ ]] && ((10#$n > last)); then
            last=$((10#$n))
        fi
    done
    CURRDIR=$(printf 'r_%03d' $((last + 1)))
fi
DEST=$ODIR/$CURRDIR

echo "PROJECT=$PROJNO"
echo "GENOME=$GENOME"
echo "MAPQ=$MAPQ"
echo "VERSION=$VERSION ($RUN_VERSION)"
echo "ODIR=$ODIR"
echo "CURRDIR=$CURRDIR"
if [[ -e "$DEST" ]]; then
    echo "NOTE: $DEST exists; files are added to it"
fi
echo

##
# project.yaml first: a missing or unreadable README.txt stops the
# delivery before anything is copied.

Rscript "$SDIR/bin/readme2yaml.R" "$PROJNO" "$ODIR" "$CURRDIR" "$GENOME"
echo

#
# The database takes pi and invest from README.txt; the delivery folder
# is usually named after the same two users.
#
YAML_PI=$(sed -n 's/^pi: //p' project.yaml)
YAML_INVEST=$(sed -n 's/^invest: //p' project.yaml)
if [[ "$ODIR" != "$DROOT/$YAML_PI/$YAML_INVEST/$PROJ" ]]; then
    warn "$ODIR is not $DROOT/<pi>/<invest>/$PROJ for pi=$YAML_PI invest=$YAML_INVEST of README.txt"
fi

##
# Copy

asbot mkdir -p "$DEST/atacSeq/docs"
asbot rsync -aP atacSeq "$DEST"
asbot rsync -aP "${DOCS[@]}" "$DEST/atacSeq/docs/"
if [[ ${#DIFF_FILES[@]} -gt 0 ]]; then
    asbot mkdir -p "$DEST/atacSeq/diff"
    asbot rsync -aP "${DIFF_FILES[@]}" "$DEST/atacSeq/diff/"
fi

##
# Email

EMAIL=deliveryEmail_${PROJNO}_${CURRDIR}_$(date +%y%m%d).txt

#
# A value as the replacement of sed s/// (git tags may contain / and &)
#
sedValue() {
    local v=${1//\\/\\\\}
    v=${v//\//\\/}
    printf '%s' "${v//&/\\&}"
}

fillTemplate() {
    sed -e "s/{PROJNO}/$(sedValue "$PROJNO")/g" \
        -e "s/{RUNDIR}/$(sedValue "$CURRDIR")/g" \
        -e "s/{VERSION}/$(sedValue "$VERSION")/g" \
        "$SDIR/$1"
}

echo
echo "========================================================================="
echo
{
    fillTemplate tmpldeliveryEmail_01.txt
    if [[ ${#DIFF_FILES[@]} -gt 0 ]]; then
        fillTemplate tmpldeliveryEmail_02.txt
    fi
    fillTemplate tmpldeliveryEmail_99.txt
} | tee "$EMAIL"
echo "========================================================================="
echo

##
# Permissions database

PERMS=("$BIC_DELIVERY/authorization_db/init_impact_project_permissions.py" -p project.yaml)

#
# Run with set +eu: the module scripts reference unset variables, and
# a failure here must say what was already done.
#
registerRun() {
    set +eu
    if [ "$(type -t module)" != "function" ]; then
        . /etc/profile.d/modules.sh || return 1
    fi
    module purge || return 1
    module load python/3.8.0 || return 1
    local python38
    python38=$(which python3.8) || return 1
    module load py-python-ldap/3.4.2 || return 1
    "$python38" "${PERMS[@]}"
}

if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRYRUN: module load python/3.8.0 py-python-ldap/3.4.2"
    echo "DRYRUN: python3.8 ${PERMS[*]@Q}"
else
    rc=0
    registerRun || rc=$?
    set -eu
    if [[ $rc -ne 0 ]]; then
        echo
        echo "    FATAL ERROR: the files are delivered to $DEST and the email"
        echo "    is in $EMAIL, but registering the run with the permissions"
        echo "    database failed (exit $rc). Fix the cause and run, here:"
        echo
        echo "        module load python/3.8.0 py-python-ldap/3.4.2"
        echo "        python3.8 ${PERMS[*]@Q}"
        echo
        exit 1
    fi
fi

echo
if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRY RUN: would deliver to $DEST; nothing was written there or to the database"
else
    echo "Delivered to $DEST"
fi
echo "Email: $EMAIL"
for w in "${WARNINGS[@]}"; do
    echo "    WARNING: $w"
done
echo
