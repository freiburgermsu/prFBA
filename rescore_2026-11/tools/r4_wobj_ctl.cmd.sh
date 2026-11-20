# R4.3 wobj control at HEAD on the PRE-REPAIR inputs (cells ffc9d3a7, _cache_panel_ec 185ca60f, 10-14 prFBA
# predictions), same pins as the main run (MSD e8446e3e, ModelSEEDpy c2e31a0, templates c12run): rebuild the
# population and re-run the 3 amplicons of the 11-15 subset control; compare with the committed (moved-aside) files.
set -o pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
PR=/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/region_validation/shared_reference
OUT=/home/freiburger/Documents/prFBA_rescore_2026-11/control/r4/wobj_ctl_prerepair
cd $SR
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
/home/freiburger/Documents/py_venv/bin/python -W ignore $SR/scripts/weighted_objective_test.py \
  --outdir $OUT --expect-modelseedpy c2e31a0 \
  --msd /home/freiburger/Documents/ModelSEEDDatabase-e8446e3e \
  --cells $PR/data/cells.csv.gz --panel-ec $PR/data/_cache_panel_ec.json --prfba $PR/data/prfba \
  --workers 3 --maxtasksperchild 2 --only 1381.4__V1-V3,1582364.3__V1-V3,1175486.3__V1-V2
rc=$?; echo "EXIT=$rc"; exit $rc
