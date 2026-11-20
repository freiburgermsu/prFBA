#!/bin/bash
# R2.2 W8 on FM: old vs new MAG cell tables of every arm (old = start-tag worktree for tracked tables,
# $W/pre_repair for the git-ignored ones; new = the installed re-score). attrib_diff --table mag.
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; PY=/home/freiburger/Documents/py_venv/bin/python
O=$W/control/w8/mag; mkdir -p $O; S=$W/start; P=$W/pre_repair
cat_cells() { local out=$1 d=$2; head -1 "$d/b000/mag_cells.csv" > $out; for b in $(ls -d $d/b* | sort); do tail -n+2 $b/mag_cells.csv >> $out; done; }
cat_cells $O/old_n500.csv  $S/mag_validation/out_batches_n500;   cat_cells $O/new_n500.csv  $R/mag_validation/out_batches_n500
cat_cells $O/old_lowc.csv  $P/mag_validation/lowc_out_batches;   cat_cells $O/new_lowc.csv  $R/mag_validation/lowc_out_batches
cat_cells $O/old_km.csv    $S/mag_validation_km/out_batches;     cat_cells $O/new_km.csv    $R/mag_validation_km/out_batches
cat_cells $O/old_rnd.csv   $P/mag_validation_rnd/out_batches;    cat_cells $O/new_rnd.csv   $R/mag_validation_rnd/out_batches
cat_cells $O/old_main_batches.csv $P/mag_validation/out_batches; cat_cells $O/new_main_batches.csv $R/mag_validation/out_batches
cmp -s $O/old_main_batches.csv $S/mag_validation/out/mag_cells_all.csv && echo "old main batches == start-tag mag_cells_all.csv" || echo "WARNING old main batches != start-tag mag_cells_all.csv"
cmp -s $O/new_main_batches.csv $R/mag_validation/out/mag_cells_all.csv && echo "new main batches == installed mag_cells_all.csv" || echo "WARNING new main batches != installed mag_cells_all.csv"
run() { local name=$1 old=$2 new=$3; echo "== $name"; $PY $W/tools/attrib_diff.py --table mag --old $old --new $new --out $O/attrib_$name.json; }
run main  $S/mag_validation/out/mag_cells_all.csv   $R/mag_validation/out/mag_cells_all.csv
run tau99 $S/mag_validation/out/mag_cells_tau99.csv $R/mag_validation/out/mag_cells_tau99.csv
run n500  $O/old_n500.csv $O/new_n500.csv
run lowc  $O/old_lowc.csv $O/new_lowc.csv
run km    $O/old_km.csv   $O/new_km.csv
run rnd   $O/old_rnd.csv  $O/new_rnd.csv
run smoke $S/mag_validation/smoke/out/mag_cells.csv $R/mag_validation/smoke/out/mag_cells.csv
