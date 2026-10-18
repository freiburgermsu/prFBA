#!/bin/bash
set -uo pipefail
MV=~/Documents/prFBA/mag_validation; MVK=~/Documents/prFBA/mag_validation_rnd
PY=~/Documents/py_venv/bin/python
TRUTH=~/Documents/prFBA/region_validation_random/data/benchmark_truth.json   # built with the manifest-taxid build_truth.py
PG=~/Documents/prFBA/random_panel/pgfam_panel_random.json
one() {
  B=$(basename "$1" .tsv)
  [ -s "$MVK/out_batches/$B/mag_cells.csv" ] && { echo "[$B] cached"; return 0; }
  "$PY" "$MV"/score_mag.py --hits "$MVK/hits/$B.tsv" --manifest "$MVK/manifests/$B.json" \
      --outdir "$MVK/out_batches/$B" --truth "$TRUTH" --pgfams "$PG" \
      > "$MVK/logs/score_$B.log" 2>&1 && echo "[$B] ok" || echo "[$B] FAILED"
}
export -f one; export MV MVK PY TRUTH PG
mkdir -p "$MVK/out_batches"
ls -1 "$MVK"/hits/*.tsv | xargs -P 5 -I{} bash -c 'one "$@"' _ {}
echo "=== RANDOM MAG SCORING DONE ==="
