#!/usr/bin/env bash
# R3.1: W3 control on the OLD (09-17) GFS, then the real refresh of data/ and data_family/, then (only
# after both caches are regenerated) move P2C kfold_audit/panel_ec_sets.json aside (stale-copy #14/#15 order).
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
W=/home/freiburger/Documents/prFBA_rescore_2026-11; R=/home/freiburger/Documents/prFBA
SR=$R/region_validation/shared_reference; S=$SR/scripts; C=$W/control/r3
OLDGFS=$W/pre_repair/region_validation/picrust2_comparison/data/genome_function_sets.json
NEWGFS=$R/region_validation/picrust2_comparison/data/genome_function_sets.json
echo "== control (old GFS $(sha256sum $OLDGFS | cut -c1-12)) $(date -Is)"
for d in data data_family; do
  $PY $S/refresh_panel_ec_cache.py --functions $OLDGFS --data-dir $SR/$d \
      --out $C/w3ctl_$d.json --expect-equal $SR/$d/_cache_panel_ec.json --report $C/w3ctl_$d.refresh.json
done
echo "== back up the old caches $(date -Is)"
for d in data data_family; do
  mkdir -p $W/pre_repair/region_validation/shared_reference/$d
  cp -p $SR/$d/_cache_panel_ec.json $W/pre_repair/region_validation/shared_reference/$d/_cache_panel_ec.json
done
echo "== real refresh (new GFS $(sha256sum $NEWGFS | cut -c1-12)) $(date -Is)"
for d in data data_family; do
  $PY $S/refresh_panel_ec_cache.py --functions $NEWGFS --data-dir $SR/$d --report $C/w3_real_$d.refresh.json
done
echo "== move kfold_audit/panel_ec_sets.json aside $(date -Is)"
K=region_validation/picrust2_comparison/data/kfold_audit
mkdir -p $W/pre_repair/$K
mv -n $R/$K/panel_ec_sets.json $W/pre_repair/$K/panel_ec_sets.json
ls -la $W/pre_repair/$K/ ; test ! -e $R/$K/panel_ec_sets.json
echo "== done $(date -Is)"
