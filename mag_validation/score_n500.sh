#!/bin/bash
set -uo pipefail
MV=~/Documents/prFBA/mag_validation
PY=~/Documents/py_venv/bin/python
one() {
  B=$(basename "$1" .tsv)
  [ -s "$MV/out_batches_n500/$B/mag_cells.csv" ] && { echo "[$B] cached"; return 0; }
  "$PY" "$MV"/score_mag.py --hits "$MV/hits_n500/$B.tsv" --manifest "$MV/manifests_n500/$B.json" \
      --outdir "$MV/out_batches_n500/$B" \
      > "$MV/logs/score_n500_$B.log" 2>&1 && echo "[$B] ok" || echo "[$B] FAILED"
}
export -f one; export MV PY
ls -1 "$MV"/hits_n500/*.tsv | xargs -P 5 -I{} bash -c 'one "$@"' _ {}
echo "=== ORIG n500 SCORING DONE ==="
