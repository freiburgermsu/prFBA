set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; CTL=/home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl; C=/home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20j_copy.json
st() { stat -c '[stat] %n size=%s mtime=%.6Y' "$@"; }
cp -p $W/control/genome_gene_families.pre_repair_d20june.json $C
sha256sum $C
st $C $R/bvbrc_cache/genome_gene_families.json
echo "== C-L4 (SNAP20J copy) -> C-D20 control selection $(date -Is)"
$PY $R/select_references.py --hits $CTL/data/hits/asv_top20_alignment_hits.json    --out $CTL/data/selection.json.tmp --gene-provider auto --gene-cache $C --select-mode legacy
mv $CTL/data/selection.json.tmp $CTL/data/selection.json
st $C $R/bvbrc_cache/genome_gene_families.json
sha256sum $C
$PY $W/tools/cmp_selection.py --old $R/region_validation/data/selection.json --new $CTL/data/selection.json    --d20 --detail-max 100000 --out $W/control/d20_code_delta.json --keys-out $W/control/d20_code_delta.keys.json
$PY $W/tools/d20_delta_detail.py $R/region_validation/data/selection.json $CTL/data/selection.json $W/control/d20_code_delta_substantive.json
echo "== C-L4J done $(date -Is)"
