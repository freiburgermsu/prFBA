#!/bin/bash
# MAG -> skani reference hits for the EmilyKin 291 dereplicated MAGs, vs
#   1) GTDB r226 (pre-sketched skani DB from the skani authors)
#   2) the BV-BRC sludge_genomes set (sketched locally; hits feed the
#      select_references.py -> build_synthetic_genomes.py prFBA downstream)
# Idempotent: sketches only if missing, searches with --resume.
set -euo pipefail
cd "$(dirname "$0")"

PY=~/Documents/py_venv/bin/python
MAGS=~/Documents/EmilyKin/dereplicated_genomes
OUT=~/Documents/EmilyKin/mag_skani_hits
GTDBTK=~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv
GTDB_DB=skani_db/skani_gtdb_r226-v0.3      # consolidated sketches.db/index.db/markers.bin
THREADS=${THREADS:-60}

# BV-BRC DB (2,161 sludge_genomes; one file is a cached BV-BRC 406 error -> skipped)
[ -d skani_db/bvbrc_sludge ] || $PY mag_skani_hits.py sketch \
    --fasta-dir ~/Documents/sludge_genomes --out skani_db/bvbrc_sludge --threads "$THREADS"

$PY mag_skani_hits.py search --mags "$MAGS" \
    --db skani_db/bvbrc_sludge --db-name bvbrc_sludge --db-kind bvbrc \
    --taxonomy ~/Documents/sludge/gtdb_taxonomy.json \
    --gtdbtk "$GTDBTK" --outdir "$OUT" --topk 20 --threads "$THREADS" --resume

$PY mag_skani_hits.py search --mags "$MAGS" \
    --db "$GTDB_DB" --db-name gtdb_r226 --db-kind gtdb \
    --taxonomy skani_db/gtdb_r226/bac120_taxonomy_r226.tsv.gz \
    --taxonomy skani_db/gtdb_r226/ar53_taxonomy_r226.tsv.gz \
    --gtdbtk "$GTDBTK" --outdir "$OUT" --topk 20 --threads "$THREADS" --resume

echo "outputs in $OUT"
