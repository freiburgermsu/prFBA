#!/bin/bash
# R4.3 objective_wrapper_equivalence.py on TF (runbook §6): ModelSEEDpy-prfba (09856f5) first on PYTHONPATH; Stage O's
# biochemistry is re-read from data/objective_tuning/meta.json (args.msd = ModelSEEDDatabase-5389554a).
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
cd $SR
[ -e data/objective_tuning/wrapper_equivalence.json ] && { echo "old wrapper_equivalence.json in place"; exit 2; }
[ "$(readlink /home/freiburger/Documents/ModelSEEDDatabase)" = ModelSEEDDatabase-5389554a ] || { echo "MSD symlink not at 5389554a"; exit 2; }
git -C /home/freiburger/Documents/ModelSEEDpy-prfba rev-parse HEAD; git -C /home/freiburger/Documents/ModelSEEDpy-prfba status --porcelain | head
t0=$(date +%s)
PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba python -W ignore scripts/objective_wrapper_equivalence.py
rc=$?; echo "wrapper rc=$rc $(( $(date +%s) - t0 )) s"; exit $rc
