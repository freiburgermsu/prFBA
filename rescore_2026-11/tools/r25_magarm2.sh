#!/usr/bin/env bash
# R2.5 step 3, second half (mag_direct_feed.py completed in rescore-r25-magarm before its OOM).
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
cd /home/freiburger/Documents/prFBA/region_validation/picrust2_comparison
s() { echo "[mag] $* ($(date -Is))"; "$@"; }
s $PY scripts/mag_prfba_ec.py
s $PY scripts/mag_score.py
echo "[mag] done $(date -Is)"
