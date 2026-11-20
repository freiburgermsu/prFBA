#!/bin/bash
# R4.3 calibrate_pfunc on TF (inside prfba-ns via tf_unit + unit_run_tf.sh). Committed run's args (scores.json
# meta.args, b9fa5c7) except --cache-dir (fresh, TF $W). Old outputs must already be in TF $W/pre_repair.
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference; W=/home/freiburger/Documents/prFBA_rescore_2026-11
cd $SR
for f in data/calibration/scores.json data/calibration/reliability.csv data/calibration/assertions.json; do
  [ -e $f ] && { echo "old output still in place: $f"; exit 2; }; done
C=$W/calib_cache_r4; [ -e $C ] && { echo "cache dir $C exists: refusing to reuse"; exit 2; }; mkdir -p $C
git -C /home/freiburger/Documents/prFBA rev-parse HEAD
sha256sum data/cells.csv.gz data/_cache_panel_ec.json ../picrust2_comparison/data/{function_maps,pgfam_products,prfba_selection,genome_function_sets}.json ../picrust2_comparison/data/comparison_cells.csv.gz
t0=$(date +%s)
python -u $SR/scripts/calibrate_pfunc.py --modes exact_tie,legacy,graduated,exact_tie_guards --densities 1.0,0.5,0.25 \
  --forced --spaces ec,rxn_ec,rxn_role,pgfam,product --bootstrap 2000 --seed 1729 --eps 0.01,0.001 \
  --outdir /home/freiburger/Documents/prFBA/region_validation/shared_reference/data/calibration/ --cache-dir $C --workers 4
rc=$?; echo "calibrate rc=$rc $(( $(date +%s) - t0 )) s"; exit $rc
