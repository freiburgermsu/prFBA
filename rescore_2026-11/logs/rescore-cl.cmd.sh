set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
st() { stat -c '[stat] %n size=%s mtime=%.6Y' "$@"; }
st $R/bvbrc_cache/genome_gene_families.json
for pair in "1 region_validation_d100_legacy" "2 region_validation_kmeans_d100_legacy" "3 region_validation_random_d100_legacy"; do
  set -- $pair; n=$1; T=$2
  echo "== C-L$n $T $(date -Is)"
  st $W/control/cl_snap_copy.json
  $PY $R/select_references.py --hits $R/$T/data/hits/asv_top20_alignment_hits.json      --out $W/control/sel_$T.json --gene-provider auto --gene-cache $W/control/cl_snap_copy.json --select-mode legacy
  st $W/control/cl_snap_copy.json $R/bvbrc_cache/genome_gene_families.json
  echo "== C-L$n compare $(date -Is)"
  $PY $W/tools/cmp_selection.py --old $R/$T/data/selection.json --new $W/control/sel_$T.json      --out $W/control/cl${n}_$T.json --keys-out $W/control/cl${n}_$T.keys.json
done
sha256sum $W/control/cl_snap_copy.json
echo "== C-L4 region_validation (d20, SNAP20 copy) -> C-D20 control selection $(date -Is)"
st /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json
$PY $R/select_references.py --hits /home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl/data/hits/asv_top20_alignment_hits.json    --out /home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl/data/selection.json.tmp --gene-provider auto --gene-cache /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json --select-mode legacy
mv /home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl/data/selection.json.tmp /home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl/data/selection.json
st /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json $R/bvbrc_cache/genome_gene_families.json
sha256sum /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json
$PY $W/tools/cmp_selection.py --old $R/region_validation/data/selection.json --new /home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl/data/selection.json    --d20 --detail-max 100000 --out $W/control/d20_code_delta.json --keys-out $W/control/d20_code_delta.keys.json
echo "== C-L chain done $(date -Is)"
