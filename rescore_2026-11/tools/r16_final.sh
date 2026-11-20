#!/usr/bin/env bash
# R1.6 cache-final checkpoint: re-stamp; if the cache moved, diff key sets (additions only expected), prove no repaired
# entry changed, and rebuild the E3 catalogue (old one copied to pre_repair first). Then FREEZE.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
$PY $W/tools/stamp_state.py --out $W/stamps/R1.6.json --keys-out $W/stamps/R1.6.cache_keys.txt    --repaired-values-out $W/stamps/R1.6.repaired_values.json --compare $W/stamps/R1.0.json | tee $W/logs/r16_stamp.out
if cmp -s $W/stamps/R1.0.repaired_values.json $W/stamps/R1.6.repaired_values.json; then echo "[PASS] the 1,006 repaired entries are unchanged"; else echo "[FAIL] a repaired entry changed"; fi
if grep -q 'CACHE CHANGED' $W/logs/r16_stamp.out; then
  echo "[info] cache moved: key-set diff"
  echo "removed: $(comm -23 <(sort $W/stamps/R1.0.cache_keys.txt) <(sort $W/stamps/R1.6.cache_keys.txt) | wc -l)"
  comm -13 <(sort $W/stamps/R1.0.cache_keys.txt) <(sort $W/stamps/R1.6.cache_keys.txt) > $W/stamps/R1.6.added_keys.txt
  echo "added: $(wc -l < $W/stamps/R1.6.added_keys.txt)"
  mkdir -p $W/pre_repair/bvbrc_cache; cp -pn $R/bvbrc_cache/e3_clade_catalogue.json $W/pre_repair/bvbrc_cache/e3_clade_catalogue.json
  $PY $R/pfunc_estimator.py build-catalogue --cache $R/bvbrc_cache/genome_gene_families.json --out $R/bvbrc_cache/e3_clade_catalogue.json
  $PY -c "import sys; sys.path.insert(0,'$R'); import pfunc_estimator as P; print('catalogue_is_current after rebuild:', P.catalogue_meta_is_current(P.catalogue_meta()))"
else
  echo "[PASS] cache unchanged since R1.0: no catalogue rebuild needed (catalogue_is_current stays true)"
fi
echo "== R1.6 FREEZE $(date -Is)"
