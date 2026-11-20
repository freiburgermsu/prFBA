#!/usr/bin/env bash
# New-side reps dumps for the six d100 trees (W8), on a copy of the post-R1.3 cache; each dump's scored
# tables must be byte-identical to the live P7 output. Then attrib_diff per tree (old = pre_repair copy).
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
C=$W/control/cache_post_r13.json
if [ ! -s $C ]; then cp -p $R/bvbrc_cache/genome_gene_families.json $C; fi
sha256sum $C $R/bvbrc_cache/genome_gene_families.json
one() { T=$1; L=${2:-}
  echo "== reps new $T $L $(date -Is)"
  $PY $W/tools/ctl_run.py --common $R/$T/scripts --cache $C --record $W/control/w8/reps_${T}_new.ctl.json --      $W/tools/dump_genecap_reps.py --out $W/control/w8/reps_${T}_new.pkl $L --keep-tables $W/control/w8/tables_${T}_new
  echo "rc=$?"
  for f in gene_capture_cells.csv gene_capture_by_region.json gene_capture_long.csv; do
    if cmp -s $W/control/w8/tables_${T}_new/$f $R/$T/data/$f; then echo "[repro-new] $T $f IDENTICAL"; else echo "[repro-new] $T $f DIFFERS"; fi
  done
  OLDP=$W/control/w8/reps_${T}_old.pkl; [ $T = region_validation_d100 ] && OLDP=$W/control/w8/reps_d100_snap.pkl
  $PY $W/tools/attrib_diff.py --table genecap --old $W/pre_repair/$T/data/gene_capture_cells.csv --new $R/$T/data/gene_capture_cells.csv      --reps-pickle $OLDP --reps-pickle $W/control/w8/reps_${T}_new.pkl --out $W/control/w8/attrib_${T}.json
}
one region_validation_d100
one region_validation_d100_legacy --legacy
one region_validation_kmeans_d100
one region_validation_kmeans_d100_legacy --legacy
one region_validation_random_d100
one region_validation_random_d100_legacy --legacy
echo "== reps new d100 done $(date -Is)"
