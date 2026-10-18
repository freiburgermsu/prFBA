#!/bin/bash
# One batch: synthesize -> search -> discard. Built to run K-way concurrently under
# xargs -P, because a single skani search saturates only ~5 cores regardless of -t,
# so concurrency ACROSS batches is what actually fills a 64-core box.
set -euo pipefail
B="$1"; DB="$2"; HITS="$3"; MANI="$4"; WORK="$5"; THREADS="${6:-6}"; SYNTHJOBS="${7:-4}"
# SCRIPTS = where synth_mags.py lives (shared). PANEL = where THIS panel's contigs,
# donors and logs live. They were one variable, which silently pointed a second panel
# at the first panel's contigs -- caught before launch, never run that way.
SCRIPTS=${SCRIPTS:-~/Documents/prFBA/mag_validation}
PANEL=${PANEL:-$SCRIPTS}
MV="$SCRIPTS"
PY=~/Documents/py_venv/bin/python
SKANI=~/Documents/prFBA/bin/skani
GRID=${GRID:-'[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0],[1.0,0.1],[0.9,0.1],[0.7,0.1],[0.5,0.1],[0.3,0.1]]'}
# Persistence depth. 100 was the original run; the depth probe (out/DEPTH_PROBE.md)
# shows selections are fully converged by 200, so 500 is sufficient with margin.
NHITS=${NHITS:-100}
TAG=$(basename "$B")
if [ -s "$HITS/$TAG.tsv" ] && [ -s "$MANI/$TAG.json" ]; then echo "[$TAG] cached"; exit 0; fi
T0=$SECONDS
rm -rf "$WORK/$TAG"; mkdir -p "$WORK/$TAG/mags"
split -d -a3 -n l/"$SYNTHJOBS" "$B" "$WORK/$TAG/chunk_" 2>/dev/null || cp "$B" "$WORK/$TAG/chunk_000"
for C in "$WORK/$TAG"/chunk_*; do
  [ -s "$C" ] || continue
  "$PY" "$SCRIPTS"/synth_mags.py --contigs "$PANEL"/contigs --donor-contigs "$PANEL"/donor_contigs \
    --donors "$PANEL"/donor_assignment.json --grid "$GRID" --ids "$C" \
    --outdir "$WORK/$TAG/mags" --manifest "$C.manifest.json" >/dev/null &
done
wait
"$PY" - "$WORK/$TAG" "$MANI/$TAG.json" <<'PYEOF'
import glob, json, sys
d, out = sys.argv[1], sys.argv[2]
rows, miss = [], []
for p in sorted(glob.glob(d + "/chunk_*.manifest.json")):
    m = json.load(open(p)); rows += m["rows"]; miss += m.get("missing", [])
json.dump(dict(n_mags=len(rows), n_missing=len(miss), missing=miss, rows=rows), open(out, "w"))
PYEOF
NM=$(ls -1 "$WORK/$TAG/mags"/*.fna 2>/dev/null | wc -l); SYN=$((SECONDS-T0)); T1=$SECONDS
ls -1 "$WORK/$TAG/mags"/*.fna > "$WORK/$TAG/q.txt"
"$SKANI" search -d "$DB" --ql "$WORK/$TAG/q.txt" -o "$HITS/$TAG.tsv.part" \
    -t "$THREADS" -n "$NHITS" --min-af 0 2>"$PANEL/logs/$TAG.skani.log"
mv "$HITS/$TAG.tsv.part" "$HITS/$TAG.tsv"
ROWS=$(( $(wc -l < "$HITS/$TAG.tsv") - 1 ))
QRET=$(cut -f2 "$HITS/$TAG.tsv" | tail -n+2 | sort -u | wc -l)
rm -rf "$WORK/$TAG"
MISS=""; [ "$QRET" -ne "$NM" ] && MISS="  !!MISSING"
echo "[$TAG] mags=$NM synth=${SYN}s search=$((SECONDS-T1))s rows=$ROWS returned=$QRET/$NM$MISS"
