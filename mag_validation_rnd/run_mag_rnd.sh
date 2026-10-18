#!/bin/bash
# Panel-as-sketch MAG arm for the UNIFORM-RANDOM panel -- a line-for-line mirror of
# mag_validation_km/run_mag_km.sh so the only thing that differs is the panel:
#   0% contamination only (the original's contaminated cells are withdrawn, so 0% is the
#   only cross-panel-comparable axis), -n 500 (probe-validated depth), batches of 500,
#   8 concurrent batches x 6 threads, synth 4-way.
# The reference sketch uses skani 0.3.2 defaults: verified to reproduce the k-means
# panel sketch (identical ANI/AF on a probe query; 108.1 vs 106.5 bytes per kbp).
set -euo pipefail
MVR=~/Documents/prFBA/mag_validation_rnd
MV=~/Documents/prFBA/mag_validation
SKANI=~/Documents/prFBA/bin/skani
export GRID='[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0]]'
export NHITS=500
ls -1 "$MVR"/contigs/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort > "$MVR"/all_ids.txt
ls -1 "$MVR"/contigs/*.fna > "$MVR"/ref_files.txt
if [ ! -s "$MVR/skani_db_panel/sketches.db" ]; then
  echo "=== sketching $(wc -l < "$MVR"/ref_files.txt) genomes (skani defaults) ==="
  "$SKANI" sketch -l "$MVR"/ref_files.txt -o "$MVR"/skani_db_panel -t 61 2> "$MVR"/logs/sketch.log
fi
echo "=== random MAG arm: $(wc -l < "$MVR"/all_ids.txt) genomes x 5 tiers, -n $NHITS ==="
rm -f "$MVR"/batches/b*; split -d -a3 -l 500 "$MVR"/all_ids.txt "$MVR"/batches/b
ls -1 "$MVR"/batches/b* | xargs -P 8 -I{} env GRID="$GRID" NHITS="$NHITS" \
  SCRIPTS="$MV" PANEL="$MVR" \
  bash "$MV"/run_batch.sh {} "$MVR"/skani_db_panel "$MVR"/hits "$MVR"/manifests "$MVR"/work 6 4
echo "=== RANDOM MAG ARM DONE: $(ls "$MVR"/hits/*.tsv | wc -l) hit files ==="
