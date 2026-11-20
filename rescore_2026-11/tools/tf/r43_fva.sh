#!/bin/bash
# R4.3 Stage O FVA variant on TF: the args recorded in data/objective_tuning_fva/meta.json (b52984f; relative paths
# as recorded) + the D11 pin --msd ModelSEEDDatabase-5389554a. ARMS / FVA_ARMS env override only if the
# selection changed (data-driven; reported).
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
MSD=/home/freiburger/Documents/ModelSEEDDatabase-5389554a
ARMS=${ARMS:-two_tier/lin/G0/X100,kn/lin/G0/X100,est/lin/G0/X1,est_shuffled/lin/G0/X1,two_tier/lin/G0/X1}
FVA_ARMS=${FVA_ARMS:-est/lin/G0/X1,kn/lin/G0/X100,two_tier/lin/G0/X1,est_shuffled/lin/G0/X1}
cd $SR
for f in data/objective_tuning_fva/per_amplicon.jsonl data/objective_tuning_fva/meta.json; do [ -e $f ] && { echo "$f exists: refusing"; exit 2; }; done
[ "$(readlink /home/freiburger/Documents/ModelSEEDDatabase)" = ModelSEEDDatabase-5389554a ] || { echo "MSD symlink not at 5389554a"; exit 2; }
echo "arms $ARMS"; echo "fva-arms $FVA_ARMS"
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
t0=$(date +%s)
python -W ignore $SR/scripts/objective_tuning.py --outdir data/objective_tuning_fva --weights data/pfunc_estimators/stage_o_weights_E3.json \
  --half confirm --arms $ARMS --fva-arms $FVA_ARMS --expect-modelseedpy c2e31a0 --msd $MSD --workers 10 --maxtasksperchild 1
rc=$?; echo "fva rc=$rc $(( $(date +%s) - t0 )) s"; exit $rc
