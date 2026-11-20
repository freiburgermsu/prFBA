#!/bin/bash
# R4.3 on TF: objective_tuning_results.py with the FVA pass (as the committed 11-03 11:27 summary), then the Stage O
# independent verification (extract, c1..c13). Needs data/objective_tuning{,_fva}/{per_amplicon.jsonl,meta.json}.
set -uo pipefail
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference; W=/home/freiburger/Documents/prFBA_rescore_2026-11
cd $SR
for f in data/objective_tuning/summary.json data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md; do [ -e $f ] && { echo "$f exists"; exit 2; }; done
[ "$(wc -l < data/objective_tuning_fva/per_amplicon.jsonl)" -ge 193 ] || { echo "FVA rows missing"; exit 2; }
t0=$(date +%s); python scripts/objective_tuning_results.py --fva-dir data/objective_tuning_fva > $W/logs/r43_otune_results.log 2>&1
rc=$?; echo "results rc=$rc $(( $(date +%s) - t0 )) s"; [ $rc = 0 ] || exit $rc
bash $W/tools/r43_otune_verify.sh
