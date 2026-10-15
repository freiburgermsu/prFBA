#!/bin/bash
# Wait for the stage-4 rebuild (the five slow bac_d100 references) to finish, verify every
# manifest exists, then resume stages 5-6. Detached so it outlives any single tool call.
cd /home/freiburger/Documents/prFBA/region_validation/shared_reference
while kill -0 540573 2>/dev/null; do sleep 60; done
echo "[chain] stage-4 driver exited $(date +%T)" >> campaign_family.log
MISSING=0
for f in 0 1 2 3 4; do
  [ -f "data_family/refs/f0${f}_bac_d100/manifest.json" ] || MISSING=$((MISSING+1))
done
if [ "$MISSING" -ne 0 ]; then
  echo "[chain] ABORT: $MISSING bac_d100 manifest(s) still missing; not resuming" >> campaign_family.log
  exit 1
fi
echo "[chain] all 30 references present; resuming stages 5-6 $(date +%T)" >> campaign_family.log
exec ~/Documents/py_venv/bin/python scripts/run_campaign.py --root data_family --rank family \
  --folds 5 --k 5 --stages 5-6 --p2-jobs 3 --p2-processes 12 --prfba-jobs 1 >> campaign_family.log 2>&1
