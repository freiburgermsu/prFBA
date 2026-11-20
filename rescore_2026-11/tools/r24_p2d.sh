#!/usr/bin/env bash
# R2.4: depth-100 as-shipped PICRUSt2 chain on the repaired cache (D3 pin e8446e3e), then the dialect test.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
P2D=$R/region_validation_d100/picrust2_comparison; P2C=$R/region_validation/picrust2_comparison
export FMAPS_ARGS="--msd /home/freiburger/Documents/ModelSEEDDatabase-e8446e3e --mst /home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0"
export SEL_WORKERS=8   # 12 workers hit a cgroup OOM at 90G (fork copy-on-write); content does not depend on it
bash $R/region_validation_d100/run_picrust2_d100.sh
echo "[r24] dialect $(date -Is)"
$PY $P2D/scripts/annotation_dialect_test.py --picrust2 $P2C/data/full_picrust2
echo "[r24] done $(date -Is)"
