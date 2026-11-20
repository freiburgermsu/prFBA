#!/usr/bin/env bash
# R2.5 checks: (1) prfba_selection on SNAP20 (full key set = the pre-repair state of 11-13) so the
# committed (09-17) -> new diff splits into cache history (committed -> SNAP20) and data (SNAP20 -> new);
# (2) the generic as-shipped checks; (3) W8 with the cache-history pairs as a documented reason.
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
SNAP20=$W/control/genome_gene_families.pre_repair_d20.json
O=$W/control/r2/d20; mkdir -p $O
SC=$R/region_validation/picrust2_comparison/scripts
D=$R/region_validation/picrust2_comparison/data; OLD=$W/pre_repair/region_validation/picrust2_comparison/data
for M in default legacy; do
  X=""; T=d20_snap20; S=prfba_selection
  if [ $M = legacy ]; then X=--legacy; T=d20_legacy_snap20; S=prfba_selection_legacy; fi
  [ -s $W/control/r2/prfba_sel_$T.json ] || $PY $W/tools/ctl_run.py --common $R/region_validation/scripts --cache $SNAP20 \
      --record $W/control/r2/prfba_sel_$T.ctl.json -- $SC/prfba_selection.py --out $W/control/r2/prfba_sel_$T.json --workers 12 $X
  $PY $W/tools/cmp_prfba_selection.py --old $OLD/$S.json --new $W/control/r2/prfba_sel_$T.json --out $O/cmp_committed_vs_$T.json \
      --cache-keys $W/stamps/R1.6.cache_keys.txt --d20 > /dev/null; echo "[chk25] committed -> $T: $($PY -c "import json; r=json.load(open('$O/cmp_committed_vs_$T.json')); print(r['n_records_differ'], r['differences_by_kind'])")"
  $PY $W/tools/cmp_prfba_selection.py --old $W/control/r2/prfba_sel_$T.json --new $D/$S.json --out $O/cmp_${T}_vs_new.json \
      --hits $R/region_validation/data/hits/asv_top20_alignment_hits.json --cache-keys $W/stamps/R1.6.cache_keys.txt --d20 > /dev/null
  echo "[chk25] $T -> new: $($PY -c "import json; r=json.load(open('$O/cmp_${T}_vs_new.json')); print({k:r.get(k) for k in ('n_records_differ','n_touched_records','n_changed_untouched','changed_subset_of_touched')})")"
done
$PY - <<PYEOF
import json, orjson
W="$W"; OLD="$OLD"
pairs=[]
for S, T, suf in (("prfba_selection","d20_snap20",""),("prfba_selection_legacy","d20_legacy_snap20","_legacy")):
    old=orjson.loads(open(f"{OLD}/{S}.json","rb").read()); ctl=orjson.loads(open(f"{W}/control/r2/prfba_sel_{T}.json","rb").read())
    for k in old:
        if k=="_meta": continue
        for rg,m in (("include","prfba_include"),("exclude_0.987","prfba_exclude")):
            if old[k][rg]!=ctl[k][rg]: pairs.append([k, m+suf])
json.dump(pairs, open("$O/cache_history_pairs.json","w")); print("[chk25] cache-history pairs (committed vs SNAP20):", len(pairs))
PYEOF
bash $W/tools/r2_aschecks.sh d20
$PY $W/tools/attrib_diff.py --table comparison --old $OLD/comparison_cells.csv.gz --new $D/comparison_cells.csv.gz \
    --selection $OLD/prfba_selection.json --selection $D/prfba_selection.json \
    --legacy-selection $OLD/prfba_selection_legacy.json --legacy-selection $D/prfba_selection_legacy.json \
    --extra-genomes newproducts=$O/gfs_changed_newproducts.json --extra-keys cache_history_0917=$O/cache_history_pairs.json \
    --d20 --out $O/w8_comparison_final.json
echo "[chk25] done $(date -Is)"
