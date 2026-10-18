#!/bin/bash
# MAG arm for the K-MEANS panel, at the probe-validated depth.
# 0% contamination only: the original panel's contaminated cells are withdrawn
# (dose-control bug, see CORRECTIONS.md), so 0% is the only cross-panel-comparable axis.
set -euo pipefail
MVK=~/Documents/prFBA/mag_validation_km
MV=~/Documents/prFBA/mag_validation
export GRID='[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0]]'
export NHITS=500
mkdir -p "$MVK"/{hits,manifests,batches,work,logs,out}
ls -1 "$MVK"/contigs/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort > "$MVK"/all_ids.txt
echo "=== k-means MAG arm: $(wc -l < "$MVK"/all_ids.txt) genomes x 5 tiers, -n $NHITS ==="
rm -f "$MVK"/batches/b*; split -d -a3 -l 500 "$MVK"/all_ids.txt "$MVK"/batches/b
ls -1 "$MVK"/batches/b* | xargs -P 8 -I{} env GRID="$GRID" NHITS="$NHITS" \
  SCRIPTS="$MV" PANEL="$MVK" \
  bash "$MV"/run_batch.sh {} "$MVK"/skani_db_panel "$MVK"/hits "$MVK"/manifests "$MVK"/work 6 4
echo "=== K-MEANS MAG ARM DONE: $(ls "$MVK"/hits/*.tsv | wc -l) hit files ==="
