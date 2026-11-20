#!/bin/bash
# R4.3 Stage O independent verification on TF: data/objective_tuning/verification/{extract, c1..c13}.py in
# VERIFICATION.md's order (cwd = the verification dir). Logs to $W/logs/r43_otune_verify/.
set -uo pipefail
V=/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/objective_tuning/verification
L=/home/freiburger/Documents/prFBA_rescore_2026-11/logs/r43_otune_verify; mkdir -p $L
cd $V; rm -f rows.pkl prior.pkl
rc_all=0
for s in extract c1_split_health c2_c3_selection c4_c5_confirm c7_mechanism c8_robust c9_est_scale c10_prior_fva c11_tables c12_selection_stability c13_volume_split; do
  t0=$(date +%s); python -u $s.py > $L/$s.log 2>&1; rc=$?; echo "== $s rc=$rc $(( $(date +%s) - t0 )) s"; [ $rc = 0 ] || rc_all=$rc
done
exit $rc_all
