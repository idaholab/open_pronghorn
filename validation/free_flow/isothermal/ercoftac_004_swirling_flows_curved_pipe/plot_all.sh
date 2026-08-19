#!/bin/bash
# Runs every plot_bev00-*.py / plot_beh00-*.py script in this directory.
# Just loops and calls python3 on each - failures (e.g. a missing exp file)
# are reported but don't stop the rest from running.

for f in plot_bev00-*.py plot_beh00-*.py; do
    echo "=== $f ==="
    python3 "$f" || echo "  -> FAILED: $f"
done
