#!/usr/bin/env bash
# Post-run checks for an as-shipped chain (R2.4 d100 / R2.5 d20).   usage: r2_aschecks.sh d100|d20
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
TAG=$1
if [ "$TAG" = d100 ]; then T=region_validation_d100; X=""; else T=region_validation; X=--d20; fi
P=$T/picrust2_comparison; D=$R/$P/data; OLD=$W/pre_repair/$P/data; O=$W/control/r2/$TAG
mkdir -p $O
echo "[chk] $TAG start $(date -Is)"
cmp $D/bvbrc_source_genomes.json $OLD/bvbrc_source_genomes.json 2>/dev/null && echo "[chk] bvbrc_source_genomes.json unchanged" || \
  { git -C $R diff --quiet -- $P/data/bvbrc_source_genomes.json 2>/dev/null; echo "[chk] bvbrc_source_genomes.json: cmp vs pre_repair copy differs or no copy"; }
git -C $R diff --exit-code --stat -- $P/data/domain_routing.csv && echo "[chk] domain_routing.csv unchanged vs HEAD" || echo "[chk] domain_routing.csv CHANGED"
for S in prfba_selection prfba_selection_legacy; do
  $PY -c "
import orjson; m=orjson.loads(open('$D/$S.json','rb').read())['_meta']
print('[chk] $S meta', {k:m[k] for k in ('n_amplicons','n_distinct_reps','n_reps_missing_from_pgfam_cache','n_sources_missing_from_pgfam_cache','git_commit')})"
  $PY $W/tools/cmp_prfba_selection.py --old $OLD/$S.json --new $D/$S.json --out $O/seldiff_$S.json \
      --hits $R/$T/data/hits/asv_top20_alignment_hits.json --cache-keys $W/stamps/R1.6.cache_keys.txt $X > $O/seldiff_$S.txt
  $PY -c "
import json; r=json.load(open('$O/seldiff_$S.json')); print('[chk] $S old->new', {k:r.get(k) for k in ('n_records_differ','differences_by_kind','n_touched_records','n_changed_untouched','changed_subset_of_touched')})"
done
$PY $W/tools/gfs_diff.py --old-gfs $OLD/genome_function_sets.json --new-gfs $D/genome_function_sets.json \
    --old-maps $OLD/function_maps.json --new-maps $D/function_maps.json \
    --old-products $W/pre_repair/$P/data/pgfam_products.json --cache $R/bvbrc_cache/genome_gene_families.json \
    --out $O/gfs_diff.json --genomes-out $O/gfs_changed_newproducts.json $X > $O/gfs_diff.txt 2>&1; echo "[chk] gfs_diff rc=$?"; cat $O/gfs_diff.txt | head -80
$PY $W/tools/attrib_diff.py --table comparison --old $OLD/comparison_cells.csv.gz --new $D/comparison_cells.csv.gz \
    --selection $OLD/prfba_selection.json --selection $D/prfba_selection.json \
    --legacy-selection $OLD/prfba_selection_legacy.json --legacy-selection $D/prfba_selection_legacy.json \
    --extra-genomes newproducts=$O/gfs_changed_newproducts.json $X --out $O/w8_comparison.json
echo "[chk] $TAG done $(date -Is)"
