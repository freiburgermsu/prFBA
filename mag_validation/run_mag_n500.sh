#!/bin/bash
# Re-run the ORIGINAL panel's MAG hit generation at the matched depth.
#
# Why: the cross-panel comparison must hold persistence depth constant, and the
# depth probe (out/DEPTH_PROBE.md) showed -n 100 leaves ~0.09% of queries abstaining
# that a deeper search resolves. Selections are fully converged by depth 200; 500 is
# used for margin and to match the k-means arm.
#
# Writes to SIBLING directories. It must never touch hits/ or out_batches/, which are
# the committed basis for RESULTS_v2 and stay the -n 100 record.
set -euo pipefail
MV=~/Documents/prFBA/mag_validation
export GRID='[[1.0,0],[0.9,0],[0.7,0],[0.5,0],[0.3,0]]'
export NHITS=500
for d in hits_n500 manifests_n500 batches_n500 work_n500 out_batches_n500 logs; do mkdir -p "$MV/$d"; done
test -d "$MV/hits" && echo "guard: existing -n 100 hits/ present and will NOT be written ($(ls "$MV"/hits/*.tsv | wc -l) files)"
ls -1 "$MV"/contigs/*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort > "$MV"/all_ids_n500.txt
echo "=== original panel, $(wc -l < "$MV"/all_ids_n500.txt) genomes x 5 tiers, -n $NHITS ==="
rm -f "$MV"/batches_n500/b*; split -d -a3 -l 500 "$MV"/all_ids_n500.txt "$MV"/batches_n500/b
ls -1 "$MV"/batches_n500/b* | xargs -P 8 -I{} env GRID="$GRID" NHITS="$NHITS" \
  SCRIPTS="$MV" PANEL="$MV" \
  bash "$MV"/run_batch.sh {} "$MV"/skani_db_panel "$MV"/hits_n500 "$MV"/manifests_n500 "$MV"/work_n500 6 4
echo "=== ORIGINAL -n 500 DONE: $(ls "$MV"/hits_n500/*.tsv | wc -l) hit files ==="
