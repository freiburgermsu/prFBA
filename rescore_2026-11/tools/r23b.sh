#!/usr/bin/env bash
# R2.3 part B: MAG panel arm / pairings at depth 100 and the depth-20 (D2) pairings and cross-panel files.
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; O=$W/control/r2/r23
cd $R
s() { echo "[r23b] $* ($(date -Is))"; "$@" || echo "[r23b] FAILED rc=$?: $*"; }
s $PY three_panel/compare_three_panels_magpanel.py
PYTHONHASHSEED=0 s $PY two_panel/compare_mag_panels.py
s $PY mag_validation/compare_arms.py --mag-cells mag_validation/out/mag_cells_tau99.csv --amp region_validation_d100/data/gene_capture_cells.csv --out mag_validation/out/paired_tau99_d100.json
# depth 20 (D2): pairings and cross-panel files with their d20 defaults
s $PY mag_validation/compare_arms.py --mag-cells mag_validation/out/mag_cells_all.csv --out mag_validation/out/paired_comparison.json
s $PY mag_validation/crossover.py mag_validation/out/mag_cells_all.csv mag_validation/lowc_out_batches/b*/mag_cells.csv
s $PY mag_validation/compare_arms.py --mag-cells mag_validation/out/mag_cells_tau99.csv --out mag_validation/out/paired_tau99.json
s $PY three_panel/compare_three_panels.py
s $PY two_panel/compare_panels.py
s $PY two_panel/mechanism_depth.py
s $PY mag_validation/plot_results.py
s $PY mag_validation/plot_crossover.py
echo "[r23b] done $(date -Is)"
