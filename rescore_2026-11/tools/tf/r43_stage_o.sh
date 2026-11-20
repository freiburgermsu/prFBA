#!/bin/bash
# R4.3 Stage O on TF (prfba-ns via tf_unit + unit_run_tf.sh). = launch_objective_tuning.sh with
# OTUNE_WEIGHTS=stage_o_weights_E3.json, plus the D11 pin --msd ModelSEEDDatabase-5389554a (explicit, so meta.json
# records the commit; the default path's symlink also points there). Validate-only first, then the full run.
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference; W=/home/freiburger/Documents/prFBA_rescore_2026-11
MSD=/home/freiburger/Documents/ModelSEEDDatabase-5389554a
OUT=$SR/data/objective_tuning; WEIGHTS=$SR/data/pfunc_estimators/stage_o_weights_E3.json
cd $SR
for f in $OUT/per_amplicon.jsonl $OUT/meta.json; do [ -e $f ] && { echo "$f exists: refusing"; exit 2; }; done
[ "$(readlink /home/freiburger/Documents/ModelSEEDDatabase)" = ModelSEEDDatabase-5389554a ] || { echo "MSD symlink not at 5389554a"; exit 2; }
[ -s $WEIGHTS ] || { echo "weights missing"; exit 2; }
git -C /home/freiburger/Documents/prFBA rev-parse HEAD
sha256sum data/weighted_objective/population.json data/weighted_objective/per_amplicon.jsonl $WEIGHTS data/_cache_panel_ec.json data/cells.csv.gz
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
python -W ignore $SR/scripts/objective_tuning.py --outdir $OUT --weights $WEIGHTS --expect-modelseedpy c2e31a0 --msd $MSD --validate-only
rc=$?; echo "validate-only rc=$rc"; [ $rc = 0 ] || exit $rc
t0=$(date +%s)
python -W ignore $SR/scripts/objective_tuning.py --outdir $OUT --weights $WEIGHTS --expect-modelseedpy c2e31a0 --msd $MSD \
  --workers 10 --maxtasksperchild 1
rc=$?; echo "objective_tuning rc=$rc $(( $(date +%s) - t0 )) s"; exit $rc
