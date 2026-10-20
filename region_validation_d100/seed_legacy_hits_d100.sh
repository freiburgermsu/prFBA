#!/usr/bin/env bash
# Seed each *_d100_legacy tree with a COPY (never a link) of its union tree's raw hits.
# Usage: seed_legacy_hits_d100.sh [region_validation_d100 region_validation_kmeans_d100 ...]
# (default: every union tree whose alignment is done). Idempotent: an identical copy is kept.
set -euo pipefail
R=/home/freiburger/Documents/prFBA
TREES=("$@"); [ ${#TREES[@]} -eq 0 ] && TREES=(region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100)
for U in "${TREES[@]}"; do
  L="${U}_legacy"; SRC="$R/$U/data/hits/raw"; DST="$R/$L/data/hits/raw"
  if [ ! -s "$SRC/align.done" ]; then echo "[seed] $U not aligned yet -- skipped"; continue; fi
  mkdir -p "$DST"
  for f in asv_top20_alignment_hits.json asv_alignment_summary.csv method_selection.json align.log; do
    if [ -f "$DST/$f" ] && [ ! -L "$DST/$f" ] && cmp -s "$SRC/$f" "$DST/$f"; then continue; fi
    rm -f "$DST/align.done" "$DST/$f"            # sentinel goes first, so a partial copy never looks done
    cp -p --reflink=never "$SRC/$f" "$DST/$f.tmp" && mv "$DST/$f.tmp" "$DST/$f"
  done
  cp -p "$SRC/align.done" "$DST/align.done"      # sentinel last
  for f in asv_top20_alignment_hits.json method_selection.json align.done; do
    cmp -s "$SRC/$f" "$DST/$f" || { echo "[seed] FAIL: $DST/$f differs" >&2; exit 1; }
    [ "$(stat -c %h "$DST/$f")" = 1 ] || { echo "[seed] FAIL: $DST/$f is hard-linked" >&2; exit 1; }
  done
  echo "[seed] $L <- $U (copies verified)"
done
