#!/bin/bash
set -uo pipefail
MV=~/Documents/prFBA/mag_validation; PY=~/Documents/py_venv/bin/python
TAU="$1"; OUT="$MV/out_batches_tau$1"
one() {
  B=$(basename "$1" .tsv)
  [ -s "$OUT/$B/mag_cells.csv" ] && { echo "[$B] cached"; return 0; }
  "$PY" "$MV"/score_mag.py --hits "$MV/hits/$B.tsv" --manifest "$MV/manifests/$B.json" \
      --outdir "$OUT/$B" --tau "$TAU" > "$MV/logs/tau${TAU}_$B.log" 2>&1 \
    && echo "[$B] ok" || echo "[$B] FAILED"
}
export -f one; export MV PY OUT TAU
mkdir -p "$OUT"
ls -1 "$MV"/hits/*.tsv | xargs -P 6 -I{} bash -c 'one "$@"' _ {}
echo "=== TAU $TAU SCORING DONE ==="
