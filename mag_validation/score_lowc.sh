#!/bin/bash
set -uo pipefail
MV=~/Documents/prFBA/mag_validation
PY=~/Documents/py_venv/bin/python
score_one() {
  B=$(basename "$1" .tsv)
  [ -s "$MV/lowc_out_batches/$B/mag_cells.csv" ] && { echo "[$B] cached"; return 0; }
  "$PY" "$MV"/score_mag.py --hits "$MV/lowc_hits/$B.tsv" --manifest "$MV/lowc_manifests/$B.json" \
      --outdir "$MV/lowc_out_batches/$B" > "$MV/logs/lowc_score_$B.log" 2>&1 \
    && echo "[$B] ok $(wc -l < "$MV/lowc_out_batches/$B/mag_cells.csv") cells" || echo "[$B] FAILED"
}
export -f score_one; export MV PY
mkdir -p "$MV/lowc_out_batches"
ls -1 "$MV"/lowc_hits/*.tsv | xargs -P 6 -I{} bash -c 'score_one "$@"' _ {}
echo "=== LOWC SCORING DONE ==="
