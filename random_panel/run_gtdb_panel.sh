#!/bin/bash
# GTDB dense-reference MAG arm for an ALTERNATIVE panel (k-means or random), built to
# match mag_validation/run_gtdb.sh exactly on the regime the three panels share:
#   * same synth_mags.py, same SEED, same fragmenting
#   * same skani binary, same GTDB r226 sketch, same -n 100 --min-af 0
#   * completeness grid {100,90,70,50,30}% at contamination 0 ONLY -- the only regime
#     the original's combine_gtdb.py reports, and the only one the contamination
#     claims survived the synthesizer audit on
# Contamination 0 needs no donor genome (synth_mags.py only requires one when a grid
# point has contamination > 0), so donors is an empty map.
# Usage: run_gtdb_panel.sh <MAG_DIR> <CONTIG_DIR> <SUB_IDS>
set -euo pipefail
MV=$1; CONTIGS=$2; SUB=$3
SKANI=~/Documents/prFBA/bin/skani
PY=~/Documents/py_venv/bin/python
SYN=~/Documents/prFBA/mag_validation/synth_mags.py
GTDB=~/Documents/prFBA/skani_db/skani_gtdb_r226-v0.3
GRID='[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0]]'
BATCH=${BATCH:-250}; SYNTHJOBS=${SYNTHJOBS:-24}; THREADS=${THREADS:-61}
mkdir -p "$MV"/{gtdb_hits,gtdb_manifests,gtdb_batches,logs,gtdb_work,donor_contigs_none}
echo '{}' > "$MV/donors_none.json"
cp "$SUB" "$MV/gtdb_sub_ids.txt"
rm -f "$MV"/gtdb_batches/b*
split -d -a3 -l "$BATCH" "$MV/gtdb_sub_ids.txt" "$MV/gtdb_batches/b"
NB=$(ls "$MV"/gtdb_batches/b* | wc -l)
echo "=== GTDB arm $MV: $(wc -l < "$MV/gtdb_sub_ids.txt") genomes -> $NB batches, grid $GRID ==="
for B in "$MV"/gtdb_batches/b*; do
  TAG=$(basename "$B")
  if [ -s "$MV/gtdb_hits/$TAG.tsv" ] && [ -s "$MV/gtdb_manifests/$TAG.json" ]; then
    echo "[$TAG] cached, skipping"; continue; fi
  T0=$SECONDS
  rm -rf "$MV/gtdb_work/$TAG"; mkdir -p "$MV/gtdb_work/$TAG/mags"
  split -d -a3 -n l/"$SYNTHJOBS" "$B" "$MV/gtdb_work/$TAG/chunk_" 2>/dev/null || cp "$B" "$MV/gtdb_work/$TAG/chunk_000"
  for C in "$MV/gtdb_work/$TAG"/chunk_*; do
    [ -s "$C" ] || continue
    "$PY" "$SYN" --contigs "$CONTIGS" --donor-contigs "$MV/donor_contigs_none" \
      --donors "$MV/donors_none.json" --grid "$GRID" --ids "$C" \
      --outdir "$MV/gtdb_work/$TAG/mags" --manifest "$C.manifest.json" >/dev/null &
  done
  wait
  "$PY" - "$MV/gtdb_work/$TAG" "$MV/gtdb_manifests/$TAG.json" <<'PYEOF'
import glob, json, sys
d, out = sys.argv[1], sys.argv[2]
rows, miss = [], []
for p in sorted(glob.glob(d + "/chunk_*.manifest.json")):
    m = json.load(open(p)); rows += m["rows"]; miss += m.get("missing", [])
json.dump(dict(n_mags=len(rows), n_missing=len(miss), missing=miss, rows=rows), open(out, "w"))
PYEOF
  NM=$(ls -1 "$MV/gtdb_work/$TAG/mags"/*.fna 2>/dev/null | wc -l)
  SYNS=$((SECONDS-T0)); T1=$SECONDS
  ls -1 "$MV/gtdb_work/$TAG/mags"/*.fna > "$MV/gtdb_work/$TAG/q.txt"
  "$SKANI" search -d "$GTDB" --ql "$MV/gtdb_work/$TAG/q.txt" \
      -o "$MV/gtdb_hits/$TAG.tsv.part" -t "$THREADS" -n 100 --min-af 0 2>"$MV/logs/gtdb.$TAG.skani.log"
  mv "$MV/gtdb_hits/$TAG.tsv.part" "$MV/gtdb_hits/$TAG.tsv"
  ROWS=$(( $(wc -l < "$MV/gtdb_hits/$TAG.tsv") - 1 ))
  QRET=$(cut -f2 "$MV/gtdb_hits/$TAG.tsv" | tail -n+2 | sort -u | wc -l)
  rm -rf "$MV/gtdb_work/$TAG"
  echo "[$TAG] mags=$NM synth=${SYNS}s search=$((SECONDS-T1))s rows=$ROWS returned=$QRET/$NM"
done
echo "=== GTDB ARM DONE ($MV): $(ls $MV/gtdb_hits/*.tsv | wc -l) hit files ==="
