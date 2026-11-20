#!/usr/bin/env bash
# R3.4 control: fold f00_d100 re-selected on PRE-REPAIR inputs in scratch outdirs (canonical dir untouched):
#  (a) "subsets": the plan's literal control: old _cache_ec_subset/_cache_pgfam_subset kept (reuse path);
#  (b) "stream":  no subsets: streams the 09-17 GFS (pre_repair copy) and the pre-repair cache SNAP (f0a84bf1),
#      i.e. the same code path the real re-select takes, on the pre-repair snapshot.
# Both read the pre-repair manifests (pre_repair/.../refs/f00_{bac,arc}_d100/manifest.json). hits/ is a COPY and
# its directory is read-only during the run (guard: an align-stamp mismatch must fail, not re-align).
set -uo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
W=/home/freiburger/Documents/prFBA_rescore_2026-11; R=/home/freiburger/Documents/prFBA
SR=$R/region_validation/shared_reference; D=$SR/data; src=$D/prfba/f00_d100
OLDGFS=$W/pre_repair/region_validation/picrust2_comparison/data/genome_function_sets.json
SNAP=$W/control/genome_gene_families.pre_repair.json
OM=$W/pre_repair/region_validation/shared_reference/data/refs
for v in subsets stream; do
  o=$W/control/r3/prfba_ctl/$v/f00_d100
  [ -d $o/hits ] && chmod u+w $o/hits
  rm -rf $o; mkdir -p $o/hits
  cp -p $src/hits/* $o/hits/
  cp -p $src/align_stamp.json $src/ref_db.json $src/ref_id_db.json $src/queries.fna $src/midas_empty.csv $o/
  [ $v = subsets ] && cp -p $src/_cache_ec_subset.json $src/_cache_pgfam_subset.json $o/
  chmod a-w $o/hits
done
for v in subsets stream; do
  o=$W/control/r3/prfba_ctl/$v/f00_d100
  /usr/bin/time -v -o $o.time.txt $PY $SR/scripts/run_prfba_fold.py --mode fold --fold 0 --density 1.0 \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --ref-manifest $OM/f00_bac_d100/manifest.json --ref-manifest $OM/f00_arc_d100/manifest.json \
    --outdir $o --select-modes legacy,exact_tie,graduated,exact_tie_guards --k2 200 --gpu-cand 400 \
    --prefilter-k 1500 --force-backend cuda --workers 30 --functions $OLDGFS --pgfam-cache $SNAP > $o.log 2>&1 &
done
wait
for v in subsets stream; do
  o=$W/control/r3/prfba_ctl/$v/f00_d100; chmod u+w $o/hits
  echo "== $v: $(grep -E '^\[(skip|done|cache|load)' $o.log | tr '\n' ' ')"
  grep -E 'Maximum resident|Elapsed' $o.time.txt
  $PY $W/tools/cmp_prfba_pred.py $src/predictions.json $o/predictions.json --out $W/control/r3/prfba_ctl/cmp_$v.json; echo "cmp rc=$?"
done
o=$W/control/r3/prfba_ctl/stream/f00_d100
$PY - $src $o <<'PYEOF'
import sys, orjson
a, b = sys.argv[1], sys.argv[2]
for f in ("_cache_ec_subset.json", "_cache_pgfam_subset.json"):
    x = orjson.loads(open(f"{a}/{f}", "rb").read()); y = orjson.loads(open(f"{b}/{f}", "rb").read())
    print(f, "content equal to the 10-14 subset:", x == y, len(x), len(y))
PYEOF
