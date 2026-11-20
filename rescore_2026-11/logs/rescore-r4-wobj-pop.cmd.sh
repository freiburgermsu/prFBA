# R4.3 weighted-objective test, stage 1 (population only), FM. D11 pin: ModelSEEDDatabase e8446e3e (worktree),
# ModelSEEDpy c2e31a0 first on PYTHONPATH (+ --expect-modelseedpy), templates = script default c12run 21a8084.
# Equivalent to launch_wobj.sh stage 1 plus the explicit --msd.
set -o pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
cd $SR
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
/home/freiburger/Documents/py_venv/bin/python -W ignore $SR/scripts/weighted_objective_test.py \
  --outdir $SR/data/weighted_objective --expect-modelseedpy c2e31a0 \
  --msd /home/freiburger/Documents/ModelSEEDDatabase-e8446e3e --population-only
rc=$?; echo "EXIT=$rc"; exit $rc
