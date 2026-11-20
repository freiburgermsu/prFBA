#!/usr/bin/env bash
# R2.5: depth-20 as-shipped PICRUSt2 chain on the repaired cache, with the D10 fix (fetch once per
# selection file) and the D3 pin (ModelSEEDDatabase ac619c6e), then ceiling, compare, render,
# dialect and the depth probe (sample -> keys must equal the committed keys -> eval on reused k500 hits).
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
cd $R/region_validation/picrust2_comparison
s() { echo "[r25] $* ($(date -Is))"; "$@"; }
s $PY scripts/prfba_selection.py --out data/prfba_selection.json --workers 12   # plan: 24; 12 here (fork copy-on-write; content independent of it)
s $PY scripts/prfba_selection.py --legacy --out data/prfba_selection_legacy.json --workers 12   # plan: 24; 12 here (fork copy-on-write; content independent of it)
s $PY scripts/fetch_bvbrc_annotations.py pgfams --selection data/prfba_selection.json
s $PY scripts/fetch_bvbrc_annotations.py pgfams --selection data/prfba_selection_legacy.json
s $PY scripts/function_maps.py --msd /home/freiburger/Documents/ModelSEEDDatabase-ac619c6e \
      --mst /home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0
s $PY scripts/genome_function_sets.py
s $PY scripts/annotation_ceiling.py
s $PY scripts/compare.py
s $PY scripts/render_comparison.py
s $PY scripts/annotation_dialect_test.py
s $PY scripts/depth_probe_sample.py
git -C $R show HEAD:region_validation/picrust2_comparison/data/depth_probe/wiped_sample.keys.json > $W/control/r2/wiped_sample.keys.committed.json
$PY -c "
import json,sys
a=json.load(open('$W/control/r2/wiped_sample.keys.committed.json')); b=json.load(open('data/depth_probe/wiped_sample.keys.json'))
same = a==b
print('[r25] depth-probe sample keys', 'IDENTICAL to committed' if same else 'DIFFER from committed')
sys.exit(0 if same else 7)"
s $PY scripts/depth_probe_eval.py
echo "[r25] done $(date -Is)"
