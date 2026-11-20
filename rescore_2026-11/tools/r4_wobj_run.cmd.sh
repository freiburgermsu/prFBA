# R4.3 weighted-objective test, stage 2 (main run), FM. = launch_wobj.sh stage 2 (--workers 8 --maxtasksperchild 2)
# plus the explicit D11 pin --msd ModelSEEDDatabase e8446e3e. Reloads population.json (spec-checked) and resumes nothing
# (per_amplicon.jsonl was moved aside).
set -o pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
cd $SR
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
/home/freiburger/Documents/py_venv/bin/python -W ignore $SR/scripts/weighted_objective_test.py \
  --outdir $SR/data/weighted_objective --expect-modelseedpy c2e31a0 \
  --msd /home/freiburger/Documents/ModelSEEDDatabase-e8446e3e --workers 8 --maxtasksperchild 2
rc=$?; echo "EXIT=$rc"; exit $rc
