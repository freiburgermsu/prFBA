#!/usr/bin/env bash
# R1.3 step 1: move every stale POST output of the six d100 trees to $W/pre_repair/<repo-relative path>
# (collect.py reads old git-ignored files from there). Checkpoint files are MOVED (so the driver
# re-runs the stage); tables that the re-run overwrites are COPIED (so the old version survives).
# Tree figures (git-ignored, renderers skip existing files) are moved so they re-render.
#   usage: r13_moveaside.sh [--dry-run]
set -euo pipefail
R=/home/freiburger/Documents/prFBA
W=/home/freiburger/Documents/prFBA_rescore_2026-11
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
LOG=$W/logs/r13_moveaside.log
say() { echo "$*" | tee -a "$LOG"; }
aside() {  # mode(mv|cp) repo-relative-path
  local mode=$1 rel=$2 src=$R/$2 dst=$W/pre_repair/$2
  [ -e "$src" ] || { say "[absent] $rel"; return 0; }
  [ -e "$dst" ] && { say "[REFUSE] $dst already exists"; exit 4; }
  if [ $DRY = 1 ]; then say "[dry $mode] $rel"; return 0; fi
  mkdir -p "$(dirname "$dst")"
  if [ "$mode" = mv ]; then mv -n "$src" "$dst"; else cp -p "$src" "$dst"; fi
  say "[$mode] $rel  $(stat -c %s "$dst") B"
}
say "== r13_moveaside $(date -Is) dry=$DRY"
for T in region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100; do
  for f in post.done gene_capture_by_region.json; do aside mv $T/data/$f; done
  for f in gene_capture_cells.csv gene_capture_long.csv; do aside cp $T/data/$f; done
  aside mv $T/figures
done
for T in region_validation_d100_legacy region_validation_kmeans_d100_legacy region_validation_random_d100_legacy; do
  for f in post.done selection.json close_exclude_rep_gap.done taxacc_region_summary.json \
           gene_capture_by_region.json mcnemar_region_pairs.csv; do aside mv $T/data/$f; done
  for f in taxacc_per_record.jsonl concordance_long.csv asv_region_concordance.csv self_recovery.csv \
           confusion_genus.csv ambiguity.csv gene_capture_cells.csv gene_capture_long.csv \
           mcnemar_regions.json close_exclude_rep_gap.log; do aside cp $T/data/$f; done
  aside mv $T/figures
done
say "== r13_moveaside done $(date -Is)"
