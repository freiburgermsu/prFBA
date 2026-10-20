#!/usr/bin/env bash
# PICRUSt2 as-shipped comparison at depth 100, ORIGINAL panel, both selection rules.
# Writes only region_validation_d100/picrust2_comparison/{data,figures}; reuses the
# depth-independent PICRUSt2 run read-only from region_validation/picrust2_comparison/data/full_picrust2.
# Needs the union AND legacy POST of the original panel (the legacy POST prefetches the
# legacy reducer's candidates; prfba_selection.py is cache-only and would read a miss as
# an empty gene set). Reads the shared PGFam cache; never writes it.
# Resumable: the two selections are skipped when present; the BV-BRC fetches resume.
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA
SC=$R/region_validation_d100/picrust2_comparison/scripts
D=$R/region_validation_d100/picrust2_comparison/data
FIG=$R/region_validation_d100/picrust2_comparison/figures
P2=$R/region_validation/picrust2_comparison/data/full_picrust2       # read-only, depth-independent
WORKERS="${SEL_WORKERS:-12}"                                         # fork pool copies ~15 GB of hits
for T in region_validation_d100 region_validation_d100_legacy; do
  [ -s "$R/$T/data/post.done" ] || { echo "[p2] $T POST not done -- run run_post_all_d100.sh first" >&2; exit 2; }
done
mkdir -p "$D" "$FIG"
run() { echo "[p2] $* ($(date -Is))"; "$@"; }
[ -s "$D/prfba_selection.json" ] || run $PY "$SC/prfba_selection.py" --out "$D/prfba_selection.json" --workers "$WORKERS"
[ -s "$D/prfba_selection_legacy.json" ] || run $PY "$SC/prfba_selection.py" --legacy --out "$D/prfba_selection_legacy.json" --workers "$WORKERS"
run $PY "$SC/fetch_bvbrc_annotations.py" genomes --outdir "$D"
run $PY "$SC/fetch_bvbrc_annotations.py" pgfams --selection "$D/prfba_selection.json" --outdir "$D"
run $PY "$SC/fetch_bvbrc_annotations.py" pgfams --selection "$D/prfba_selection_legacy.json" --outdir "$D"
run $PY "$SC/function_maps.py" --products "$D/pgfam_products.json" --out "$D/function_maps.json"
run $PY "$SC/genome_function_sets.py" --maps "$D/function_maps.json" \
    --selection "$D/prfba_selection.json,$D/prfba_selection_legacy.json" --out "$D/genome_function_sets.json"
run $PY "$SC/compare.py" --picrust2 "$P2" \
    --selection "$D/prfba_selection.json" --legacy-selection "$D/prfba_selection_legacy.json" \
    --functions "$D/genome_function_sets.json" --maps "$D/function_maps.json" \
    --genomes "$D/bvbrc_source_genomes.json" --cells "$D/comparison_cells.csv.gz" \
    --summary "$D/comparison_summary.json" --routing "$D/domain_routing.csv"
run $PY "$SC/render_comparison.py" --summary "$D/comparison_summary.json" --outdir "$FIG"
echo "[p2] done $(date -Is) -> $D/comparison_summary.json"
