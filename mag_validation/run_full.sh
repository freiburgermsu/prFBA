#!/bin/bash
# Full MAG-arm run: synthesize -> search -> discard, in batches.
#
# Batched because the full grid is ~262 GB of synthetic FASTA if kept on disk;
# each batch is searched and then deleted, holding peak disk to ~15 GB. Resumable:
# a batch whose hits TSV already exists is skipped.
set -euo pipefail
MV=~/Documents/prFBA/mag_validation
SKANI=~/Documents/prFBA/bin/skani
PY=~/Documents/py_venv/bin/python
GRID='[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0],[1.0,0.1],[0.9,0.1],[0.7,0.1],[0.5,0.1],[0.3,0.1]]'
BATCH=${BATCH:-500}          # source genomes per batch -> BATCH*10 synthetic MAGs
SYNTHJOBS=${SYNTHJOBS:-24}   # parallel synthesizer processes
THREADS=${THREADS:-61}       # ~95% of 64 cores for skani

mkdir -p "$MV"/{hits,manifests,batches,logs,work}
ls -1 "$MV"/contigs/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort > "$MV"/all_ids.txt
TOTAL=$(wc -l < "$MV"/all_ids.txt)
split -d -a3 -l "$BATCH" "$MV"/all_ids.txt "$MV"/batches/b
NB=$(ls "$MV"/batches/b* | wc -l)
echo "=== $TOTAL genomes, $BATCH per batch -> $NB batches, 10 tiers each ==="

for B in "$MV"/batches/b*; do
  TAG=$(basename "$B")
  if [ -s "$MV/hits/$TAG.tsv" ] && [ -s "$MV/manifests/$TAG.json" ]; then
    echo "[$TAG] cached, skipping"; continue
  fi
  T0=$SECONDS
  rm -rf "$MV/work/$TAG"; mkdir -p "$MV/work/$TAG/mags"
  # parallel synthesis: split the batch across SYNTHJOBS chunks
  rm -f "$MV/work/$TAG"/chunk_*
  split -d -a3 -n l/"$SYNTHJOBS" "$B" "$MV/work/$TAG/chunk_" 2>/dev/null || cp "$B" "$MV/work/$TAG/chunk_000"
  for C in "$MV/work/$TAG"/chunk_*; do
    [ -s "$C" ] || continue
    "$PY" "$MV"/synth_mags.py --contigs "$MV"/contigs --donor-contigs "$MV"/donor_contigs \
      --donors "$MV"/donor_assignment.json --grid "$GRID" --ids "$C" \
      --outdir "$MV/work/$TAG/mags" --manifest "$C.manifest.json" >/dev/null &
  done
  wait
  "$PY" - "$MV/work/$TAG" "$MV/manifests/$TAG.json" <<'PYEOF'
import glob, json, sys
d, out = sys.argv[1], sys.argv[2]
rows, miss = [], []
for p in sorted(glob.glob(d + "/chunk_*.manifest.json")):
    m = json.load(open(p)); rows += m["rows"]; miss += m.get("missing", [])
json.dump(dict(n_mags=len(rows), n_missing=len(miss), missing=miss, rows=rows), open(out, "w"))
print("  synthesized %d MAGs (%d skipped)" % (len(rows), len(miss)))
PYEOF
  NM=$(ls -1 "$MV/work/$TAG/mags"/*.fna 2>/dev/null | wc -l)
  SYN=$((SECONDS-T0)); T1=$SECONDS
  ls -1 "$MV/work/$TAG/mags"/*.fna > "$MV/work/$TAG/q.txt"
  "$SKANI" search -d "$MV"/skani_db_panel --ql "$MV/work/$TAG/q.txt" \
      -o "$MV/hits/$TAG.tsv.part" -t "$THREADS" -n 100 --min-af 0 2>"$MV/logs/$TAG.skani.log"
  mv "$MV/hits/$TAG.tsv.part" "$MV/hits/$TAG.tsv"
  ROWS=$(( $(wc -l < "$MV/hits/$TAG.tsv") - 1 ))
  QRET=$(cut -f2 "$MV/hits/$TAG.tsv" | tail -n+2 | sort -u | wc -l)
  SRCH=$((SECONDS-T1))
  rm -rf "$MV/work/$TAG"
  echo "[$TAG] mags=$NM synth=${SYN}s search=${SRCH}s rows=$ROWS queries_returned=$QRET/$NM disk_avail=$(df -h /home/freiburger | awk 'NR==2{print $4}')"
  if [ "$QRET" -ne "$NM" ]; then echo "  !! WARNING: $((NM-QRET)) queries returned no hits at all"; fi
done
echo "=== ALL BATCHES DONE: $(ls $MV/hits/*.tsv | wc -l) hit files, $(cat $MV/hits/*.tsv | wc -l) total lines ==="
