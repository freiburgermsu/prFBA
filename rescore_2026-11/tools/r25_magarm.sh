#!/usr/bin/env bash
# R2.5 step 3: P2C MAG arm (deck only) on the new d20 GFS/maps, the rebuilt panel index and the R2.2 batch tables.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
cd /home/freiburger/Documents/prFBA/region_validation/picrust2_comparison
s() { echo "[mag] $* ($(date -Is))"; "$@"; }
s $PY scripts/mag_direct_feed.py
s $PY scripts/mag_prfba_ec.py
s $PY scripts/mag_score.py
echo "[mag] done $(date -Is)"
