#!/usr/bin/env bash
# Cross-panel comparisons at depth 100, both selection rules. Needs post.done in all six
# trees. Every output is a NEW file (suffix _d100 / _d100_legacy); no depth-20 output is
# written. Read-only on the PGFam cache.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA
for T in region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100; do
  for X in "$T" "${T}_legacy"; do
    [ -s "$R/$X/data/post.done" ] || { echo "[down] $X POST not done" >&2; exit 2; }
  done
done
run() { echo "[down] $* ($(date -Is))"; "$@"; }
for S in _d100 _d100_legacy; do
  O=$R/region_validation$S; K=$R/region_validation_kmeans$S
  run $PY "$R/three_panel/compare_three_panels.py" --suffix "$S"      # -> three_panel/three_panel_comparison$S.json
  run $PY "$R/three_panel/emptied_strata.py" --suffix "$S"            # -> three_panel/emptied_strata$S.json (EMPTIED)
  run $PY "$R/two_panel/compare_panels.py" \
      --orig-truth "$O/data/benchmark_truth.json" --km-truth "$K/data/benchmark_truth.json" \
      --orig-genecap "$O/data/gene_capture_cells.csv" --km-genecap "$K/data/gene_capture_cells.csv" \
      --out "$R/two_panel/panel_comparison$S.json"
  run $PY "$R/two_panel/mechanism_depth.py" --depth 100 --original-root "$O" --kmeans-root "$K" \
      --out "$R/two_panel/mechanism_depth$S.json"
  # self-recovery diagnosis (original panel; writes that tree's data/self_recovery_*)
  run $PY "$O/scripts/diagnose_self_recovery.py"
  # MAG arm vs the 16S arm (the cited depth-20 runs used the original panel's cells).
  # --mag-cells / the crossover inputs are inferred from the cited outputs' tier sets: verify.
  run $PY "$R/mag_validation/compare_arms.py" --amp "$O/data/gene_capture_cells.csv" \
      --mag-cells "$R/mag_validation/out/mag_cells_all.csv" \
      --out "$R/mag_validation/out/paired_comparison$S.json"
  run $PY "$R/mag_validation/crossover.py" --amp "$O/data/gene_capture_cells.csv" \
      --out "$R/mag_validation/out/crossover$S.json" \
      "$R/mag_validation/out/mag_cells_all.csv" "$R"/mag_validation/lowc_out_batches/b*/mag_cells.csv
done
run $PY "$R/region_validation_d100/check_d100.py" all
echo "[down] done $(date -Is)"
