#!/bin/bash
# Build all vignettes to HTML
# Usage: ./build_vignettes.sh [vignette_name.Rmd]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ -n "$1" ]; then
    # Build specific vignette
    echo "Building $1..."
    Rscript -e "rmarkdown::render('$1', output_format = 'html_document')"
else
    # Build all vignettes
    for f in *.Rmd; do
        echo "Building $f..."
        Rscript -e "rmarkdown::render('$f', output_format = 'html_document')"
    done
fi

echo "Done! HTML files created."
