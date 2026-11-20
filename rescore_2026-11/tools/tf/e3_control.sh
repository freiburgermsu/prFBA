#!/bin/bash
# C-E3 (plan R2.1 step 1) on TF, inside prfba-ns started with PRFBA_NS_BINDS that put the pre-repair cache +
# old catalogue over bvbrc_cache/, the cc8014f reference cells over the canonical cell paths, and
# $W/control_e3/out/{pfunc_e3_gate.json,pfunc_e3_gate/} over the gate's output paths (so the clone's tracked
# outputs are never written). Today's code (HEAD) must reproduce the committed 'results' blocks.
set -uo pipefail
R=/home/freiburger/Documents/prFBA; cd $R
ls -la --time-style=full-iso bvbrc_cache/ region_validation_d100/data/gene_capture_cells*.csv
sha256sum region_validation_d100/data/gene_capture_cells.csv region_validation_d100/data/gene_capture_cells_tau0995.csv bvbrc_cache/e3_clade_catalogue.json
git rev-parse HEAD
python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 10
echo "gate rc=$?"
