#!/bin/bash
# watch_and_run.sh REF_EPOCH
# Wait until the notebook's cell-12 merge rewrites model_inputs/BV_BRC_16S.json with a
# fresh mtime (> REF_EPOCH), the shard files are gone, and the file is stable + brace-valid;
# then embed every unique 16S sequence and analyze the resulting space.
# Writes sentinels DONE / FAILED into the store and logs everything to run.log.
set -u
REF="${1:?need REF_EPOCH}"
REPO=/home/freiburger/Documents/codiffusion_bioreactor
OUT=/home/freiburger/Documents/prFBA
PY=/home/freiburger/Documents/py_venv/bin/python
TARGET="$REPO/model_inputs/BV_BRC_16S.json"
LOG="$OUT/run.log"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True   # avoid VRAM fragmentation on the 16GB card
rm -f "$OUT/DONE" "$OUT/FAILED"

ts() { date '+%Y-%m-%d %H:%M:%S'; }
say() { echo "[$(ts)] $*" >> "$LOG"; }

say "watcher armed; waiting for $TARGET mtime > $REF ($(date -d @"$REF"))"
heartbeat=0
while :; do
    if [ -f "$TARGET" ]; then
        m=$(stat -c %Y "$TARGET" 2>/dev/null || echo 0)
        if [ "$m" -gt "$REF" ]; then
            nshards=$(ls "$REPO"/Fasta16S*.json 2>/dev/null | wc -l)
            s1=$(stat -c %s "$TARGET"); sleep 60; s2=$(stat -c %s "$TARGET")
            first=$(head -c1 "$TARGET" 2>/dev/null); last=$(tail -c1 "$TARGET" 2>/dev/null)
            say "candidate complete: mtime fresh, shards=$nshards, size ${s1}->${s2}, ends '${first}..${last}'"
            if [ "$nshards" -eq 0 ] && [ "$s1" = "$s2" ] && [ "$first" = "{" ] && [ "$last" = "}" ]; then
                break
            fi
        fi
    fi
    heartbeat=$((heartbeat + 1))
    if [ $((heartbeat % 10)) -eq 0 ]; then
        ns=$(ls "$REPO"/Fasta16S*.json 2>/dev/null | wc -l)
        say "still waiting... shards=$ns, target mtime=$(stat -c %Y "$TARGET" 2>/dev/null) (ref=$REF)"
    fi
    sleep 30
done

say "DATASET COMPLETE: $(stat -c 'size=%s mtime=%y' "$TARGET")"
say "=== embedding (embed_16s.py) ==="
cd "$OUT" || { echo FAILED > "$OUT/FAILED"; exit 1; }
if ! $PY "$OUT/embed_16s.py" --batch-size 128 >> "$LOG" 2>&1; then
    say "EMBEDDING FAILED"; echo "embed failed" > "$OUT/FAILED"; exit 1
fi
say "=== analysis (analyze_space.py) ==="
if ! $PY "$OUT/analyze_space.py" >> "$LOG" 2>&1; then
    say "ANALYSIS FAILED"; echo "analyze failed" > "$OUT/FAILED"; exit 1
fi
say "ALL DONE"
echo "ok" > "$OUT/DONE"
