#!/usr/bin/env bash
# R1.3 checks (plan R1.3 step 3), all read-only. Light (< 10 GB): the selection comparisons stream.
#   1 check_d100.py post <T> for the six trees (check 5)
#   2 predicted output: original union cells/summary/long == the 11-13 repaired-cache prediction
#   3 union trees: selection.json and taxacc_per_record.jsonl sha256 == R1.0 stamp; git diff on the tree's
#     tracked data lists only the two gene-capture files
#   4 legacy trees: record-level selection comparison old (pre_repair) vs new; changed must be touched
#   5 fetch counts in the POST log (any > 0 moves the cache stamp: list for R1.6)
#   6 headline numbers old -> new for the six gene_capture_by_region.json
set -u
R=/home/freiburger/Documents/prFBA
W=/home/freiburger/Documents/prFBA_rescore_2026-11
PY=/home/freiburger/Documents/py_venv/bin/python
OUT=$W/control/r13; mkdir -p $OUT
UNION="region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100"
LEGACY="region_validation_d100_legacy region_validation_kmeans_d100_legacy region_validation_random_d100_legacy"
echo "== 1 check_d100 post"
for T in $UNION $LEGACY; do echo "-- $T"; $PY $R/region_validation_d100/check_d100.py post $T; echo "rc=$?"; done
echo "== 2 predicted output (original union vs expected/repro_0987)"
for f in gene_capture_cells.csv gene_capture_by_region.json gene_capture_long.csv; do
  if cmp -s $R/region_validation_d100/data/$f $W/expected/repro_0987/$f; then echo "[PASS] $f byte-identical to the 11-13 prediction"
  else echo "[FAIL] $f differs from the 11-13 prediction"; fi; done
echo "== 3 union trees unchanged outside gene capture"
for T in $UNION; do
  for f in selection.json taxacc_per_record.jsonl; do
    want=$(grep " $T/data/$f\$" $W/stamps/ignored_before.sha256 | cut -d' ' -f1)
    have=$(sha256sum $R/$T/data/$f | cut -d' ' -f1)
    [ "$want" = "$have" ] && echo "[PASS] $T/data/$f sha256 == R1.0" || echo "[FAIL] $T/data/$f sha256 $have != R1.0 $want"
  done
  ch=$(cd $R && GIT_OPTIONAL_LOCKS=0 git diff --name-only -- $T/data | sort | tr '\n' ' ')
  echo "[info] $T tracked data changed: $ch"
  extra=$(echo $ch | tr ' ' '\n' | grep -vE '/gene_capture_(cells\.csv|by_region\.json)$' | grep -v '^$')
  [ -z "$extra" ] && echo "[PASS] $T: only gene-capture files changed" || echo "[FAIL] $T: unexpected changes: $extra"
done
echo "== 4 legacy selections: changed records must be touched"
for T in $LEGACY; do
  $PY $W/tools/cmp_selection.py --old $W/pre_repair/$T/data/selection.json --new $R/$T/data/selection.json \
      --detail-max 3000 --out $OUT/sel_$T.json --keys-out $OUT/sel_$T.keys.json
done
echo "== 5 fetch counts in the POST log"
grep -nE 'to_fetch=|fetched [0-9]+|missing from cache|merged [0-9]+ new|nothing to fetch' $R/region_validation_d100/logs/post_all_rescore.log
echo "== 6 headline numbers"
args=()
for T in $UNION $LEGACY; do args+=(--pair $T $W/pre_repair/$T/data/gene_capture_by_region.json $R/$T/data/gene_capture_by_region.json); done
$PY $W/tools/genecap_headline.py "${args[@]}" --out $OUT/genecap_headline_d100.json
stat -c '[stat] %n size=%s mtime=%.6Y' $R/bvbrc_cache/genome_gene_families.json
echo "== r13 checks done $(date -Is)"
