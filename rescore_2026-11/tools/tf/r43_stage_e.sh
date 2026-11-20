#!/bin/bash
# R4.3 Stage E on TF (prfba-ns via tf_unit + unit_run_tf.sh): pfunc_estimators.py with the committed run's args
# (scores.json meta.args, b46598b) except --cache-dir (fresh), then the E3 weights producer (95675ac).
# Needs: the new calibration/scores.json (r43-calib) and the new weighted_objective/population.json (pushed from FM).
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference; W=/home/freiburger/Documents/prFBA_rescore_2026-11
cd $SR
for f in data/pfunc_estimators/scores.json data/pfunc_estimators/stage_o_weights.json data/pfunc_estimators/stage_o_weights_E3.json; do
  [ -e $f ] && { echo "old output still in place: $f"; exit 2; }; done
[ -s data/weighted_objective/population.json ] || { echo "population.json missing"; exit 2; }
[ -s data/calibration/scores.json ] || { echo "calibration/scores.json missing"; exit 2; }
C=$W/pfe_cache_r4; [ -e $C ] && { echo "cache dir $C exists: refusing to reuse"; exit 2; }; mkdir -p $C
git -C /home/freiburger/Documents/prFBA rev-parse HEAD
sha256sum data/cells.csv.gz data/_cache_panel_ec.json data/weighted_objective/population.json data/calibration/scores.json
t0=$(date +%s)
python -u $SR/scripts/pfunc_estimators.py --workers 12 --bootstrap 2000 --seed 1729 --densities 1.0,0.5,0.25 \
  --outdir data/pfunc_estimators --cache-dir $C
rc=$?; echo "pfunc_estimators rc=$rc $(( $(date +%s) - t0 )) s"; [ $rc = 0 ] || exit $rc
python -u $SR/scripts/stage_o_weights_e3.py
rc=$?; echo "stage_o_weights_e3 rc=$rc"; exit $rc
