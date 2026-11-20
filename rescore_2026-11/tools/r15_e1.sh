#!/usr/bin/env bash
# R1.5 (E1): genecap_tau0995 on the repaired cache and the new P7 cells. Smoke first, then full.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
aside() { local mode=$1 rel=$2; local dst=$W/pre_repair/$rel
  [ -e "$R/$rel" ] || { echo "[absent] $rel"; return 0; }
  [ -e "$dst" ] && { echo "[REFUSE] $dst exists"; exit 4; }
  mkdir -p "$(dirname "$dst")"; if [ $mode = mv ]; then mv -n "$R/$rel" "$dst"; else cp -pr "$R/$rel" "$dst"; fi
  echo "[$mode] $rel"; }
D=region_validation_d100/data
aside mv $D/genecap_tau0995/pgfam_supplement.json
aside mv $D/genecap_tau0995/smoke
aside mv $D/genecap_tau0995/repro_0987
aside cp $D/genecap_tau0995/run_meta.json
for f in gene_capture_cells_tau0995.csv gene_capture_long_tau0995.csv gene_capture_by_region_tau0995.json gene_capture_tau_compare.json; do aside cp $D/$f; done
cd $R/region_validation_d100
echo "== E1 smoke (--limit 300) $(date -Is)"
$PY -u scripts/genecap_tau0995.py --limit 300 > logs/genecap_tau0995_rescore_smoke.log 2>&1
tail -5 logs/genecap_tau0995_rescore_smoke.log
echo "== E1 full $(date -Is)"
$PY -u scripts/genecap_tau0995.py > logs/genecap_tau0995_rescore.log 2>&1
tail -8 logs/genecap_tau0995_rescore.log
ls -la data/genecap_tau0995/
echo "== E1 done $(date -Is)"
