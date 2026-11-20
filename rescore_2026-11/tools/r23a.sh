#!/usr/bin/env bash
# R2.3 part A: the depth-100 downstream driver (D1-D8, both rules; check_d100 all is expected to fail check 6
# on the d20 baseline files rewritten by R1.4), then the hash-seed proof of compare_three_panels.
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; O=$W/control/r2/r23
cd $R
bash $R/region_validation_d100/run_downstream_d100.sh; echo "[r23a] downstream driver rc=$?"
for S in _d100 _d100_legacy; do
  PYTHONHASHSEED=777 $PY three_panel/compare_three_panels.py --suffix $S --out $O/three_panel_comparison${S}.seed777.json > $O/tpc${S}_seed777.log 2>&1
  cmp three_panel/three_panel_comparison${S}.json $O/three_panel_comparison${S}.seed777.json && echo "[r23a] compare_three_panels $S: byte-identical under PYTHONHASHSEED 0 and 777" || echo "[r23a] compare_three_panels $S: DIFFERS between hash seeds"
done
echo "[r23a] done $(date -Is)"
