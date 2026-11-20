#!/usr/bin/env bash
# R1.4 (decision D2): re-select and re-score the d20 original tree (legacy) with today's code on the
# repaired cache, and re-score gene capture in the d20 k-means/random union trees. By hand (hazard H5:
# never the d20 run_validation.sh), under the d100 cache-writer flock (taken non-blocking by the caller).
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; CACHE=$R/bvbrc_cache/genome_gene_families.json
st() { stat -c '[stat] %n size=%s mtime=%.6Y' "$@"; }
aside() { local mode=$1 rel=$2; local dst=$W/pre_repair/$rel
  [ -e "$R/$rel" ] || { echo "[absent] $rel"; return 0; }
  [ -e "$dst" ] && { echo "[REFUSE] $dst exists"; exit 4; }
  mkdir -p "$(dirname "$dst")"; if [ $mode = mv ]; then mv -n "$R/$rel" "$dst"; else cp -p "$R/$rel" "$dst"; fi
  echo "[$mode] $rel"; }
[ ! -e $CACHE.tmp ] || { echo "$CACHE.tmp exists"; exit 3; }
st $CACHE
D=region_validation/data
aside mv $D/selection.json
for f in taxacc_per_record.jsonl concordance_long.csv asv_region_concordance.csv self_recovery.csv confusion_genus.csv          ambiguity.csv taxacc_region_summary.json gene_capture_cells.csv gene_capture_long.csv gene_capture_by_region.json          mcnemar_region_pairs.csv mcnemar_regions.json self_recovery_diagnosis.json self_recovery_failures.csv; do aside cp $D/$f; done
for T in region_validation_kmeans region_validation_random; do
  for f in gene_capture_cells.csv gene_capture_long.csv gene_capture_by_region.json; do aside cp $T/data/$f; done; done
echo "== d20 candidate prefetch (legacy reducer) $(date -Is)"
$PY $R/prefetch_gene_families.py --hits $R/$D/hits/asv_top20_alignment_hits.json --cache $CACHE --workers 16 --select-mode legacy
st $CACHE
echo "== d20 select_references legacy $(date -Is)"
$PY $R/select_references.py --hits $R/$D/hits/asv_top20_alignment_hits.json    --out $R/$D/selection.json.tmp --gene-provider auto --gene-cache $CACHE --select-mode legacy
mv $R/$D/selection.json.tmp $R/$D/selection.json
st $CACHE
echo "== d20 gap close $(date -Is)"
$PY $R/kmeans_panel/close_exclude_rep_gap.py --tree region_validation --legacy 2>&1 | tee -a $R/$D/close_exclude_rep_gap.rescore.log
st $CACHE
cd $R/region_validation/scripts
echo "== d20 score_taxacc --legacy --force $(date -Is)"; $PY score_taxacc.py --legacy --force
echo "== d20 score_genecap --legacy --force $(date -Is)"; $PY score_genecap.py --legacy --force
echo "== d20 mcnemar --force $(date -Is)"; $PY mcnemar_regions.py --force
echo "== d20 diagnose_self_recovery $(date -Is)"; $PY diagnose_self_recovery.py
for T in region_validation_kmeans region_validation_random; do
  echo "== $T score_genecap --force $(date -Is)"; (cd $R/$T/scripts && $PY score_genecap.py --force)
done
st $CACHE
echo "== r14 done $(date -Is)"
