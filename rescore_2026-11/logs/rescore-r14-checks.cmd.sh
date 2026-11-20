PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
for T in region_validation_kmeans region_validation_random; do for f in selection.json taxacc_per_record.jsonl; do
  [ -e $R/$T/data/$f ] || { echo "[absent] $T/data/$f"; continue; }
  want=$(grep " $T/data/$f$" $W/stamps/ignored_before.sha256 | cut -d' ' -f1); have=$(sha256sum $R/$T/data/$f | cut -d' ' -f1)
  [ "$want" = "$have" ] && echo "[PASS] $T/data/$f sha256 == R1.0" || echo "[FAIL] $T/data/$f $have != $want"; done
  echo "[info] $T tracked data changed: $(cd $R && GIT_OPTIONAL_LOCKS=0 git diff --name-only -- $T/data | tr '\n' ' ')"; done
echo "[info] region_validation tracked data changed: $(cd $R && GIT_OPTIONAL_LOCKS=0 git diff --name-only -- region_validation/data | tr '\n' ' ')"
echo "== d20 selection: control (today's code, pre-repair SNAP20J) -> new (repaired cache), --d20"
$PY $W/tools/cmp_selection.py --old $W/control/region_validation_ctl/data/selection.json --new $R/region_validation/data/selection.json --d20 --detail-max 5000 --out $W/control/r14/sel_d20_ctl_to_new.json --keys-out $W/control/r14/sel_d20_ctl_to_new.keys.json
echo "== d20 selection: committed (June) -> new, --d20"
$PY $W/tools/cmp_selection.py --old $W/pre_repair/region_validation/data/selection.json --new $R/region_validation/data/selection.json --d20 --detail-max 2000 --out $W/control/r14/sel_d20_committed_to_new.json
args=(--pair region_validation $W/pre_repair/region_validation/data/gene_capture_by_region.json $R/region_validation/data/gene_capture_by_region.json
      --pair region_validation_ctl_to_new $W/control/region_validation_ctl/data/gene_capture_by_region.json $R/region_validation/data/gene_capture_by_region.json
      --pair region_validation_kmeans $W/pre_repair/region_validation_kmeans/data/gene_capture_by_region.json $R/region_validation_kmeans/data/gene_capture_by_region.json
      --pair region_validation_random $W/pre_repair/region_validation_random/data/gene_capture_by_region.json $R/region_validation_random/data/gene_capture_by_region.json)
$PY $W/tools/genecap_headline.py "${args[@]}" --out $W/control/r14/genecap_headline_d20.json
echo "== r14 checks done"
