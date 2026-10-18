#!/bin/bash
set -uo pipefail
MV=~/Documents/prFBA/mag_validation
PY=~/Documents/py_venv/bin/python
score_one() {
  B=$(basename "$1" .tsv)
  [ -s "$MV/out_batches/$B/mag_cells.csv" ] && { echo "[$B] cached"; return 0; }
  "$PY" "$MV"/score_mag.py --hits "$MV/hits/$B.tsv" --manifest "$MV/manifests/$B.json" \
      --outdir "$MV/out_batches/$B" > "$MV/logs/score_$B.log" 2>&1 \
    && echo "[$B] ok $(wc -l < "$MV/out_batches/$B/mag_cells.csv") cells" \
    || echo "[$B] FAILED"
}
export -f score_one; export MV PY
ls -1 "$MV"/hits/*.tsv | xargs -P 6 -I{} bash -c 'score_one "$@"' _ {}
echo "=== SCORING DONE ==="
