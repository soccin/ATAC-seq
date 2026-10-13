#!/bin/bash
#
# Render the client-facing docs (output.md, methods.md) to PDF and HTML
# with pandoc. Run from this directory. Needs pandoc and xelatex with the
# DejaVu fonts.
#

set -euo pipefail

cd "$(dirname "$0")"

#
# render MD TITLE [PANDOC_OPTS ...]
#
render() {
    local MD=$1
    local TITLE=$2
    local BASE=${MD%.md}
    shift 2

    pandoc $MD -o $BASE.pdf "$@" \
        --pdf-engine=xelatex \
        -V sansfont="DejaVu Sans" \
        -V monofont="DejaVu Sans Mono" \
        -V geometry:margin=.75in \
        -V fontsize=12pt \
        -V header-includes="\\usepackage{fancyhdr} \\pagestyle{fancy} \\fancyhead[L]{$TITLE} \\fancyhead[R]{\\today}" \
        --highlight-style=kate

    pandoc $MD -o $BASE.html \
        --standalone \
        --css=https://sindresorhus.com/github-markdown-css/github-markdown.css \
        --metadata title="$TITLE"
}

render output.md "ATAC-seq Output" --number-sections
render methods.md "ATAC-seq Methods"
