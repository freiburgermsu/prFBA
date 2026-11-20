#!/usr/bin/env bash
# R2.4/R2.5 control: prfba_selection on the pre-repair snapshots must reproduce the committed selections.
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
SNAP=$W/control/genome_gene_families.pre_repair.json
SNAP20J=$W/control/genome_gene_families.pre_repair_d20june.json
O=$W/control/r2
P2D=$R/region_validation_d100/picrust2_comparison; P2C=$R/region_validation/picrust2_comparison
rc_all=0
run_ctl() {  # tree snap workers tag extra
  local T=$1 S=$2 NW=$3 TAG=$4; shift 4
  local SC=$R/$T/picrust2_comparison/scripts
  echo "[selctl] $TAG start $(date -Is)"
  $PY $W/tools/ctl_run.py --common $R/$T/scripts --cache $S --record $O/prfba_sel_$TAG.ctl.json -- \
      $SC/prfba_selection.py --out $O/prfba_sel_$TAG.json --workers $NW "$@" || { echo "[selctl] $TAG RUN FAILED"; rc_all=1; return; }
  echo "[selctl] $TAG done $(date -Is)"
}
cmp_ctl() {  # tag committed
  $PY $W/tools/cmp_prfba_selection.py --old $2 --new $O/prfba_sel_$1.json --out $O/cmp_$1.json \
      --cache-keys $W/stamps/R1.6.cache_keys.txt && echo "[selctl] $1 RECORDS IDENTICAL" || { echo "[selctl] $1 RECORDS DIFFER"; rc_all=1; }
}
run_ctl region_validation_d100 $SNAP 12 d100 && cmp_ctl d100 $P2D/data/prfba_selection.json
run_ctl region_validation_d100 $SNAP 12 d100_legacy --legacy && cmp_ctl d100_legacy $P2D/data/prfba_selection_legacy.json
run_ctl region_validation $SNAP20J 24 d20_snap20j && cmp_ctl d20_snap20j $P2C/data/prfba_selection.json
run_ctl region_validation $SNAP20J 24 d20_legacy_snap20j --legacy && cmp_ctl d20_legacy_snap20j $P2C/data/prfba_selection_legacy.json
echo "[selctl] all done rc=$rc_all $(date -Is)"
exit $rc_all
