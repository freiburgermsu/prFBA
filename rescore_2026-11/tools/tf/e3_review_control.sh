#!/bin/bash
# Review-recompute control on TF (prfba-ns with the pre-repair cache + old catalogue bound over bvbrc_cache/):
# today's pfunc_e3_gate_review_recompute.py on the committed per_amplicon.csv.gz must reproduce the
# committed review_recompute.json (all but elapsed_s).
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; cd $R
ls -la --time-style=full-iso bvbrc_cache/; sha256sum region_validation_d100/data/pfunc_e3_gate/per_amplicon.csv.gz; git rev-parse HEAD
mkdir -p $W/control_e3/review
PFUNC_GATE_REVIEW_OUT=$W/control_e3/review WORKERS=8 python -u region_validation_d100/scripts/pfunc_e3_gate_review_recompute.py > $W/control_e3/review/stdout.json
echo "review rc=$?"
