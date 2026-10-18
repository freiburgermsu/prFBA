#!/bin/bash
# Dense-reference arm: the same synthetic MAGs searched against GTDB r226 (732,475
# accessions) instead of the 10,000-genome panel.
#
# Why a second arm: the panel is species-dereplicated, so after self-exclusion the
# nearest surviving relative is a genus-level one sitting near skani's ~83 ANI
# detection floor -- completeness cannot move an abstention that distance already
# decides. A dense reference has marginal-AF relatives whose survival completeness
# DOES decide, which is the frontier this arm measures.
#
# Scope: a phylum-stratified 2,000-genome subsample (20,000 queries), because GTDB
# search costs ~0.4 s/query vs ~0.03 s against the panel.
# Metrics here need NO accession mapping: ANI, aligned fraction, novelty class,
# and abstention under include / exclude(ANI>=95). Taxonomy and gene capture stay
# in the panel arm, where every reference is a BV-BRC genome with NCBI lineage and
# cached PGFams.
set -euo pipefail
MV=~/Documents/prFBA/mag_validation
SKANI=~/Documents/prFBA/bin/skani
PY=~/Documents/py_venv/bin/python
GTDB=~/Documents/prFBA/skani_db/skani_gtdb_r226-v0.3
GRID=${GRID:-'[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0],[1.0,0.1],[0.9,0.1],[0.7,0.1],[0.5,0.1],[0.3,0.1]]'}
BATCH=${BATCH:-250}
SYNTHJOBS=${SYNTHJOBS:-24}
THREADS=${THREADS:-61}

mkdir -p "$MV"/{gtdb_hits,gtdb_manifests,gtdb_batches,logs,gtdb_work}
rm -f "$MV"/gtdb_batches/b*
split -d -a3 -l "$BATCH" "$MV"/gtdb_sub_ids.txt "$MV"/gtdb_batches/b
NB=$(ls "$MV"/gtdb_batches/b* | wc -l)
echo "=== GTDB arm: $(wc -l < "$MV"/gtdb_sub_ids.txt) genomes, $BATCH per batch -> $NB batches ==="

for B in "$MV"/gtdb_batches/b*; do
  TAG=$(basename "$B")
  if [ -s "$MV/gtdb_hits/$TAG.tsv" ] && [ -s "$MV/gtdb_manifests/$TAG.json" ]; then
    echo "[$TAG] cached, skipping"; continue
  fi
  T0=$SECONDS
  rm -rf "$MV/gtdb_work/$TAG"; mkdir -p "$MV/gtdb_work/$TAG/mags"
  split -d -a3 -n l/"$SYNTHJOBS" "$B" "$MV/gtdb_work/$TAG/chunk_" 2>/dev/null || cp "$B" "$MV/gtdb_work/$TAG/chunk_000"
  for C in "$MV/gtdb_work/$TAG"/chunk_*; do
    [ -s "$C" ] || continue
    "$PY" "$MV"/synth_mags.py --contigs "$MV"/contigs --donor-contigs "$MV"/donor_contigs \
      --donors "$MV"/donor_assignment.json --grid "$GRID" --ids "$C" \
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
  SYN=$((SECONDS-T0)); T1=$SECONDS
  ls -1 "$MV/gtdb_work/$TAG/mags"/*.fna > "$MV/gtdb_work/$TAG/q.txt"
  "$SKANI" search -d "$GTDB" --ql "$MV/gtdb_work/$TAG/q.txt" \
      -o "$MV/gtdb_hits/$TAG.tsv.part" -t "$THREADS" -n 100 --min-af 0 2>"$MV/logs/gtdb.$TAG.skani.log"
  mv "$MV/gtdb_hits/$TAG.tsv.part" "$MV/gtdb_hits/$TAG.tsv"
  ROWS=$(( $(wc -l < "$MV/gtdb_hits/$TAG.tsv") - 1 ))
  QRET=$(cut -f2 "$MV/gtdb_hits/$TAG.tsv" | tail -n+2 | sort -u | wc -l)
  rm -rf "$MV/gtdb_work/$TAG"
  echo "[$TAG] mags=$NM synth=${SYN}s search=$((SECONDS-T1))s rows=$ROWS returned=$QRET/$NM"
done
echo "=== GTDB ARM DONE: $(ls $MV/gtdb_hits/*.tsv | wc -l) hit files ==="
