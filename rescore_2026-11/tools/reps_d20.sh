#!/usr/bin/env bash
# d20 W8 (after R1.4): control-side dump on the C-D20 control tree (SNAP20J copy), new-side dumps of the three
# d20 trees on a copy of the post-R1.4 cache; attribution control -> new (d20 original, --d20) and
# committed -> new (all three; d20 original also with --d20).
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; CTL=/home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl
C=$W/control/cache_post_r14.json
if [ ! -s $C ]; then cp -p $R/bvbrc_cache/genome_gene_families.json $C; fi
sha256sum $C $R/bvbrc_cache/genome_gene_families.json
cmpt() { for f in gene_capture_cells.csv gene_capture_by_region.json gene_capture_long.csv; do
    if cmp -s $1/$f $2/$f; then echo "[repro] $3 $f IDENTICAL"; else echo "[repro] $3 $f DIFFERS"; fi; done; }
echo "== reps control (C-D20 tree, SNAP20J) $(date -Is)"
$PY $W/tools/ctl_run.py --common $CTL/scripts --cache $W/control/d20ctl/snap20j_copy.json --record $W/control/w8/reps_d20ctl.ctl.json --    $W/tools/dump_genecap_reps.py --out $W/control/w8/reps_d20ctl.pkl --legacy --keep-tables $W/control/w8/tables_d20ctl
cmpt $W/control/w8/tables_d20ctl $CTL/data d20ctl
for spec in "region_validation --legacy" "region_validation_kmeans" "region_validation_random"; do
  set -- $spec; T=$1; L=${2:-}
  echo "== reps new $T $L $(date -Is)"
  $PY $W/tools/ctl_run.py --common $R/$T/scripts --cache $C --record $W/control/w8/reps_${T}_new.ctl.json --      $W/tools/dump_genecap_reps.py --out $W/control/w8/reps_${T}_new.pkl $L --keep-tables $W/control/w8/tables_${T}_new
  cmpt $W/control/w8/tables_${T}_new $R/$T/data $T
done
echo "== W8 d20 original: control -> new (--d20)"
$PY $W/tools/attrib_diff.py --table genecap --d20 --old $CTL/data/gene_capture_cells.csv --new $R/region_validation/data/gene_capture_cells.csv    --reps-pickle $W/control/w8/reps_d20ctl.pkl --reps-pickle $W/control/w8/reps_region_validation_new.pkl --out $W/control/w8/attrib_region_validation_ctl_to_new.json
echo "== W8 d20 original: committed -> new (--d20; informative)"
$PY $W/tools/attrib_diff.py --table genecap --d20 --old $W/pre_repair/region_validation/data/gene_capture_cells.csv --new $R/region_validation/data/gene_capture_cells.csv    --reps-pickle $W/control/w8/reps_region_validation_old.pkl --reps-pickle $W/control/w8/reps_d20ctl.pkl --reps-pickle $W/control/w8/reps_region_validation_new.pkl    --out $W/control/w8/attrib_region_validation_committed_to_new.json
for T in region_validation_kmeans region_validation_random; do
  echo "== W8 $T committed -> new"
  $PY $W/tools/attrib_diff.py --table genecap --old $W/pre_repair/$T/data/gene_capture_cells.csv --new $R/$T/data/gene_capture_cells.csv      --reps-pickle $W/control/w8/reps_${T}_old.pkl --reps-pickle $W/control/w8/reps_${T}_new.pkl --out $W/control/w8/attrib_${T}.json
done
echo "== reps d20 done $(date -Is)"
