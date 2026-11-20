#!/bin/bash
# R2.1 on TF (prfba-ns, canonical repaired inputs): smoke twice (determinism), the gate, the review recompute.
# The old outputs must already be in $W/pre_repair (moved aside before this unit).
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; cd $R
D=region_validation_d100/data
for f in $D/pfunc_e3_gate.json $D/pfunc_e3_gate/per_amplicon.csv.gz $D/pfunc_e3_gate/run_meta.json $D/pfunc_e3_gate/review_recompute.json; do
  [ -e $f ] && { echo "old output still in place: $f"; exit 2; }; done
[ -e $D/pfunc_e3_gate/smoke ] && { echo "smoke dir exists"; exit 2; }
ls -la --time-style=full-iso bvbrc_cache/; sha256sum $D/gene_capture_cells.csv $D/gene_capture_cells_tau0995.csv bvbrc_cache/e3_clade_catalogue.json; git rev-parse HEAD
for i in 1 2; do
  t0=$(date +%s)
  python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --limit 500 --workers 8 --verify-q > $W/logs/e3_smoke$i.log 2>&1
  rc=$?; echo "smoke$i rc=$rc $(( $(date +%s) - t0 )) s"; rm -rf $W/e3_smoke$i; mv $D/pfunc_e3_gate/smoke $W/e3_smoke$i
done
t0=$(date +%s)
python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 10
rc=$?; echo "gate rc=$rc $(( $(date +%s) - t0 )) s"
[ $rc = 0 ] || exit $rc
t0=$(date +%s); mkdir -p $W/e3_review
PFUNC_GATE_REVIEW_OUT=$W/e3_review WORKERS=8 python -u region_validation_d100/scripts/pfunc_e3_gate_review_recompute.py > $W/e3_review/stdout.json
rc=$?; echo "review rc=$rc $(( $(date +%s) - t0 )) s"
[ $rc = 0 ] && cp $W/e3_review/indep_recompute.json $D/pfunc_e3_gate/review_recompute.json
exit $rc
