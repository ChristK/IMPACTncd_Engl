#!/bin/bash
# Watch for changes to Rmd files and automatically rebuild them
# Requires: inotify-tools (sudo dnf install inotify-tools)
#
# Usage: ./watch_vignettes.sh
# Stop with Ctrl+C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Watching for changes to .Rmd files in $SCRIPT_DIR"
echo "Press Ctrl+C to stop"
echo ""

# Check if inotifywait is available
if ! command -v inotifywait &> /dev/null; then
    echo "Error: inotifywait not found. Install with:"
    echo "  sudo dnf install inotify-tools   # Rocky/RHEL/Fedora"
    echo "  sudo apt install inotify-tools   # Debian/Ubuntu"
    exit 1
fi

# Watch for file modifications
inotifywait -m -e close_write --format '%w%f' . | while read FILE; do
    if [[ "$FILE" == *.Rmd ]]; then
        echo ""
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Changed: $FILE"
        echo "Rendering..."
        
        Rscript -e "rmarkdown::render('$FILE', output_format = 'html_document', quiet = TRUE)" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ Successfully rendered $(basename "${FILE%.Rmd}.html")"
        else
            echo "✗ Error rendering $FILE"
        fi
    fi
done
